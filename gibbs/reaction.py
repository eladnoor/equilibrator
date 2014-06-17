# -*- coding: utf-8 -*-

import logging
import numpy
import urllib
import os

from gibbs import conditions
from gibbs import constants
from gibbs import models

relpath = os.path.dirname(os.path.realpath(__file__))
cc_preprocess_fname = os.path.join(relpath, '../data/cc_preprocess.npz')
cc_preprocess = numpy.load(cc_preprocess_fname)

class ReactantFormulaMissingError(Exception):
    
    def __init__(self, c):
        self.compound = c
        
    def __str__(self):
        return ("Cannot test reaction balancing because the reactant %s does not "
                "have a chemical formula" % self.compound.kegg_id)

class CompoundWithCoeff(object):
    """A compound with a stoichiometric coefficient."""
    
    def __init__(self, coeff, compound, phase, name=None):
        """Construct a CompoundWithCoeff object.
        
        Args:
            coeff: the coefficient.
            compound: the Compound object.
            name: a string of the compound name.
            concentration: the concentration (molar).
        """
        self.compound = compound
        self.coeff = coeff
        self.phase = phase
        self._name = name
        self.transformed_energy = None
    
    def Clone(self):
        other = CompoundWithCoeff(self.coeff, self.compound,
                                  self.phase, self.name)
        return other
    
    def AbsoluteCoefficient(self):
        return numpy.abs(self.coeff)
    
    @staticmethod
    def FromReactant(reactant):
        return CompoundWithCoeff(reactant.coeff, reactant.compound,
                                 conditions.StandardAqueousPhase(),
                                 name=reactant.compound.FirstName())
    
    @staticmethod
    def FromId(coeff, kegg_id, phase=None, name=None):
        try:
            compound = models.Compound.objects.get(kegg_id=kegg_id)
        except Exception:
            return None
        if name is None:
            name = name or compound.FirstName()
        if phase is None:
            def_phase_name = compound.GetDefaultPhaseName()
            phase = conditions._BaseConditions._GeneratePhase(def_phase_name, 1e-3)
        
        return CompoundWithCoeff(coeff, compound, phase, name)

    @staticmethod
    def FromDict(d):
        kegg_id = d['kegg_id']
        try:
            compound = models.Compound.objects.get(kegg_id=kegg_id)
        except Exception:
            return None

        coeff = d.get('coeff', 1)
        name = d.get('name', compound.FirstName())
        phase_name = d.get('phase', None)
        if phase_name is None:
            phase_name = compound.GetDefaultPhaseName()
        conc = d.get('conc', 1e-3)

        phase = conditions._BaseConditions._GeneratePhase(phase_name, conc)

        logging.debug('Compound = %s, name = %s, coeff = %d, %s' %
                      (kegg_id, name, coeff, phase))
        return CompoundWithCoeff(coeff, compound, phase, name)
        
    def Minus(self):
        """Returns a new CompoundWithCoeff with coeff = -self.coeff."""
        return CompoundWithCoeff(-self.coeff, self.compound, self.phase,
                                 self.name)
    
    def __str__(self):
        name = self.name or str(self.compound)
        return '%d %s' % (self.coeff, name)
    
    def ToJson(self, include_species=True):
        d = {'coeff': self.coeff,
             'KEGG_ID': self.compound.kegg_id,
             'phase': str(self.phase),
             'name': str(self.compound.first_name),
             'source_used': None}
        
        if self.compound._species_group is not None:
            d['source_used'] = str(self.compound._species_group.formation_energy_source)
        
        if include_species:
            d['species'] = self.compound.SpeciesJson()
            
        return d
        
    def _GetUrlParams(self):
        return ['reactantsId=%s' % self.kegg_id,
                'reactantsCoeff=%d' % self.coeff,
                'reactantsName=%s' % self.name,
                'reactantsPhase=%s' % self.phase.PhaseName(),
                'reactantsConcentration=%s' % self.phase.Value()]
    
    def GetCompoundList(self):
        sg = self.compound._species_group
        species = sg.GetPhaseSpecies(self.phase.PhaseName())
        logging.debug('KEGG ID = %s' % self.kegg_id)
        logging.debug('# species in phase %s = %d' % (self.phase.PhaseName(), len(species)))
        
        l = []
        for s in species:
            l.append({'nh': int(s.number_of_hydrogens),
                      'charge': int(s.net_charge),
                      'nmg': int(s.number_of_mgs),
                      'dgzero': float(s.formation_energy)})
        logging.debug('compound_list = %s' % str(l))
        return l
    
    def GetName(self):
        """Gives a string name for this compound."""
        if self.compound.preferred_name:
            return self.compound.preferred_name
        if self._name:
            return self._name
        return str(self.compound.FirstName())
    
    def __eq__(self, other):
        """
            Check equality with another CompoundWithCoeff.
        
            Args:
                other: a second CompoundWithCoeff (or like object).
        """
        if self.coeff != other.coeff:
            return False
        
        if self.GetKeggID() != other.GetKeggID():
            return False
            
        if self.phase.Name() != other.phase.Name():
            return False
        
        return True

    def DeltaG0Prime(self, aq_params):
        dg0_prime = self.compound.DeltaG0Prime(aq_params,
                                               self.phase.PhaseName())
        return self.coeff * dg0_prime
        
    def GetPhaseName(self):
        return self.phase.PhaseName()
        
    def GetPossiblePhaseNames(self):
        return self.compound.GetPossiblePhaseNames()
        
    def HasMultiplePhases(self):
        return len(self.GetPossiblePhaseNames()) > 1

    def GetPhaseSubscript(self):
        # The compound C00288 represents the group of all carbonate species
        # and CO2 and is called CO2(total). It is technically aqueous, but
        # also contains a gas component, so the phase is not well defined
        # and hence the (total) subscript. In this case, we don't need to
        # add another subscript so we return an empty string
        if self.name == 'CO2(total)':
            return ''
        return self.phase.Subscript()
    
    def GetPhaseValueString(self):
        return '%.8g' % self.phase.HumanValueAndUnits()[0]

    def GetPhasePrefactor(self):
        return self.phase.HumanValueAndUnits()[1]
        
    def GetPhaseUnits(self):
        return self.phase.HumanValueAndUnits()[2]

    def GetPhaseIsConstant(self):
        return self.phase.IsConstant()
        
    def GetPhaseHumanString(self):
        value, prefactor, units = self.phase.HumanValueAndUnits_letters()
        return '%g %s' % (value, units)

    def GetKeggID(self):
        return self.compound.kegg_id
    
    name = property(GetName)
    abs_coeff = property(AbsoluteCoefficient)
    subscript = property(GetPhaseSubscript)
    phase_name = property(GetPhaseName)
    human_string = property(GetPhaseHumanString)
    phase_value_string = property(GetPhaseValueString)
    phase_prefactor = property(GetPhasePrefactor)
    phase_units = property(GetPhaseUnits)
    is_constant = property(GetPhaseIsConstant)
    possible_phases = property(GetPossiblePhaseNames)
    has_multiple_phases = property(HasMultiplePhases)
    kegg_id = property(GetKeggID)

class Reaction(object):
    """A reaction."""
    
    def __init__(self, reactants=None, aq_params=None):
        """Construction.
        
        Args:
            substrates: a list of CompoundWithCoeff objects.
            products: a list of CompoundWithCoeff objects.
        """
        self.reactants = reactants or []
        self.aq_params = aq_params or conditions.AqueousParams()

        self._FilterProtons()
        self._Dedup()

        self._kegg_id = None

        # used only as cache, no need to copy while cloning        
        self._catalyzing_enzymes = None
        
    def SetAqueousParams(self, aq_params):
        self._dg0_prime = None
        self._aq_params = aq_params
        self._SetCompoundPriorities()
    
    def GetAqueousParams(self):
        return self._aq_params
        
    aq_params = property(GetAqueousParams, SetAqueousParams)

    def Clone(self):
        """
            Make a copy of the reaction
        """
        logging.debug('Cloning reaction...')
        other = Reaction([r.Clone() for r in self.reactants], 
                         self.aq_params.Clone())
        other._kegg_id = self._kegg_id
        return other

    def _GetMaxCommonPriority(self):
        max_priority = self.aq_params.max_priority or 1
        
        # The chosen priority will be the highest number which is common
        # to all reactants.
        priorities = [c.compound.GetSpeciesGroupPriorities() for c in self.reactants]
        priorities = filter(lambda l: len(l) > 0, priorities)
        
        # Someone is missing data!
        if priorities == []:
            return 0
        else:
            return min([max(l) for l in priorities] + [max_priority])

    def _SetCompoundPriorities(self):
        """
            Returns a set of (int, SpeciesGroup) tuples for the reaction.
        """
        priority = self._GetMaxCommonPriority()
        for c in self.reactants:
            c.compound.SetSpeciesGroupPriority(priority)
         
    def SwapSides(self):
        """Swap the sides of this reaction."""
        for c in self.reactants:
            c.coeff = -c.coeff
        
    def GetSubstrates(self):
        s = [c for c in self.reactants if c.coeff < 0]
        return s
        
    def GetProducts(self):
        p = [c for c in self.reactants if c.coeff > 0]
        return p
    
    def __str__(self):
        """
            Simple text reaction representation.
        """
        return self.GetQueryString()

    def GetSparseRepresentation(self):
        sparse = {}
        for r in self.reactants:
            sparse[r.kegg_id] = r.coeff
        return sparse

    def GetHashableReactionString(self):
        return models.StoredReaction.HashableReactionString(self.GetSparseRepresentation())
    
    def GetHash(self):
        return models.StoredReaction.HashReaction(self.GetSparseRepresentation())
    
    def _GetAllStoredReactions(self):
        """
            Find all stored reactions matching this compound (using the hash)
        """
        logging.debug('looking for stored reactions matching this one')
        my_hash = self.GetHash()
        my_string = self.GetHashableReactionString()

        matching_stored_reactions = models.StoredReaction.objects.select_related(
            ).filter(reaction_hash=my_hash)
        logging.debug('my hash = %s (%d matches)' %
                      (my_hash, len(matching_stored_reactions)))
        
        matching_stored_reactions = [m for m in matching_stored_reactions
            if m.GetHashableReactionString() == my_string]
        logging.debug('my hashable string = %s (%d matches)' %
                      (my_string, len(matching_stored_reactions)))
        
        return matching_stored_reactions
    
    def _GetCatalyzingEnzymes(self):
        """
            Get all the enzymes catalyzing this reaction.
        """
        if self._catalyzing_enzymes is None:
            logging.debug('looking for enzymes catalyzing this reaction')
            self._catalyzing_enzymes = []
            for stored_reaction in self._GetAllStoredReactions():
                enzymes = stored_reaction.enzyme_set.all()
                self._catalyzing_enzymes.extend(enzymes)
        return self._catalyzing_enzymes
    
    def ToJson(self):
        """
            Return this reaction as a JSON-compatible object.
        """
        cdicts = [c.ToJson(include_species=False) for c in self.reactants]
        enzdicts = [e.ToJson() for e in self._GetCatalyzingEnzymes()]
        d = {'reaction_string': str(self), 
             'reactants': cdicts,
             'enzymes': enzdicts,
             'chemically_balanced': self.is_balanced,
             'redox_balanced': self.is_electron_balanced,
             'dgzero': None,
             'dgzero_prime': None,
             'keq_prime': None,
             'KEGG_ID': self._kegg_id}
        
        if self.dg0_prime is not None:
            d['dgzero_prime'] = {
                'value': round(self.dg0_prime, 1),
                'pH': self.aq_params.pH,
                'ionic_strength': self.aq_params.ionic_strength}
        if self.k_eq_prime is not None:
            d['keq_prime'] = {
                'value': self.k_eq_prime,
                'pH': self.aq_params.pH,
                'ionic_strength': self.aq_params.ionic_strength}
        
        return d

    @staticmethod
    def FromForm(form):
        """
            Build a reaction object from a ReactionForm.
            
            Args:
                form: a ReactionForm object.
            
            Returns:
                A Reaction object or None if there's an error.
        """
        max_priority = form.cleaned_max_priority
        n_react = len(form.cleaned_reactantsCoeff)
        coeffs = list(form.cleaned_reactantsCoeff)
        kegg_ids = list(form.cleaned_reactantsId)
        names = list(form.cleaned_reactantsName)
        phases = list(form.cleaned_reactantsPhase)
        concentrations = list(form.cleaned_reactantsConcentration)

        compound_list = []
        for i in xrange(n_react):
            d = {'coeff': coeffs[i], 'kegg_id': kegg_ids[i], 'name': names[i]}
            if phases != []:
                d['phase'] = phases[i]
            if concentrations != []:
                logging.debug("concentrations = %s" % str(concentrations))
                d['conc'] = concentrations[i]
            compound_list.append(d)

        # Return the built reaction object.
        return Reaction.FromIds(compound_list)
    
    @staticmethod
    def FromIds(compound_list):
        """Build a reaction object from lists of IDs.
        
        Args:
            compound_list: an iterable of (coeff, kegg_id, name) of reactants.
            concentration_profile: a ConcentrationProfile object.
            
        Returns:
            A properly set-up Reaction object or None if there's an error.
        """        
        reactants = map(CompoundWithCoeff.FromDict, compound_list)
        return Reaction(reactants)
        
    @staticmethod
    def _GetCollectionAtomDiff(collection):
        """Get the net atom counts from the collection.
        
        Args:
            collection: an iterable of CompoundWithCoeff instances.
        """
        atom_diff = {}
        for compound_w_coeff in collection:
            c = compound_w_coeff.compound
            coeff = compound_w_coeff.coeff

            atom_bag = c.GetAtomBag()
            if not atom_bag:
                logging.warning('Failed to fetch atom bag for %s', c.formula)
                raise ReactantFormulaMissingError(c)
            
            for atomic_number, atom_count in atom_bag.iteritems():
                new_diff = atom_diff.get(atomic_number, 0) - coeff * atom_count
                atom_diff[atomic_number] = new_diff
        
        return atom_diff
    
    @staticmethod
    def _GetCollectionElectronDiff(collection):
        """Get the net electron count from the collection.
        
        Args:
            collection: an iterable of CompoundWithCoeff instances.
        """
        electron_diff = 0
        for compound_w_coeff in collection:
            c = compound_w_coeff.compound
            coeff = compound_w_coeff.coeff
            
            electrons = c.num_electrons
            if electrons == None:
                logging.warning('Compound %s has unknown electron count', 
                                c.kegg_id)
                return 0
            
            electron_diff += coeff * electrons
        
        return electron_diff
    
    def _GetAtomDiff(self):
        """Returns the net atom counts from this reaction."""
        return self._GetCollectionAtomDiff(self.reactants)
    
    def _GetElectronDiff(self):
        """Returns the net electron count from this reaction."""
        return self._GetCollectionElectronDiff(self.reactants)
    
    @staticmethod
    def _IsBalanced(atom_diff):
        """Checks if the per-atom diffs represent a balanced collection.
        
        Args:
            atom_diff: a dictionary mapping atomic numbers to counts.
            
        Returns:
            True if balanced.
        """
        # Always ignore hydrogens, ala Alberty.
        atom_diff.pop('H', 0)
        
        # This can happen for reactions that only involve hydrogen
        if len(atom_diff) == 0:
            return True
        
        return max([abs(x) for x in atom_diff.values()]) < 0.01
    
    def _GetUrlParams(self, query=None):
        """
            Get the URL params for this reaction.
        """
        params = sum([c_w_c._GetUrlParams() for c_w_c in self.reactants], [])
        params.extend(self.aq_params._GetUrlParams())
                
        if query is not None:
            for arrow in constants.POSSIBLE_REACTION_ARROWS:
                tmp_query = query.replace(arrow, '=>')
            params.append('query=%s' % urllib.quote(tmp_query))
        
        return params
    
    def GetHyperlink(self, query=None):
        params = self._GetUrlParams(query)
        return '/reaction?%s' % '&'.join(params)
    
    def GetBalanceWithWaterLink(self, query=None):
        """Returns a link to balance this reaction with water."""
        other = self.Clone()
        other.TryBalanceWithWater()
        return other.GetHyperlink(query)

    def GetBalanceElectronsLink(self, query=None):
        """
            Returns a link to the same reaction,
            balanced by extra NAD+/NADH pairs.
        """
        other = self.Clone()
        other.BalanceElectrons()
        return other.GetHyperlink(query)

    def GetReplaceCO2Link(self, query=None):
        """
            Returns a link to the same reaction, but with HCO3-
            instead of CO2.
        """
        other = self.Clone()
        co2_id = 'C00011'
        bic_id = 'C00288'
        other._ReplaceCompound(co2_id, bic_id)        
        other.TryBalanceWithWater()
        link = other.GetHyperlink(query)
        return link

    def GetTemplateData(self, query=None):
        template_data = {'reaction': self,
                         'query': query}
        try:
            balance_with_water_link = self.GetBalanceWithWaterLink(query)
            balance_electrons_link = self.GetBalanceElectronsLink(query)
            replace_co2_link = self.GetReplaceCO2Link(query)
            template_data.update({
                'balance_with_water_link': balance_with_water_link,
                'balance_electrons_link': balance_electrons_link,
                'replace_co2_link': replace_co2_link})
        except ReactantFormulaMissingError:
            pass
        template_data.update(self.aq_params.GetTemplateData())
        return template_data
    
    def GetNewPriorityLink(self, max_priority):
        params = self._GetUrlParams()
        params.append('max_priority=%d' % max_priority)
        return '/reaction?%s' % '&'.join(params)
    
    def GetPhGraphLink(self):
        params = self._GetUrlParams()
        params.append('vary_ph=1')
        return '/graph_reaction?%s' % '&'.join(params)

    def GetPMgGraphLink(self):
        params = self._GetUrlParams()
        params.append('vary_pmg=1')
        return '/graph_reaction?%s' % '&'.join(params)
    
    def GetIonicStrengthGraphLink(self):
        params = self._GetUrlParams()
        params.append('vary_is=1')
        return '/graph_reaction?%s' % '&'.join(params)
    
    @staticmethod
    def _GetReactionSideString(side):
        """Write a reaction side as a string."""
        sdata = []
        for c_w_coeff in side:
            if c_w_coeff.coeff == 1:
                sdata.append(c_w_coeff.GetName())
            else:
                sdata.append('%d %s' % (c_w_coeff.coeff,
                                        c_w_coeff.GetName()))
        return ' + '.join(sdata)
    
    def GetQueryString(self):
        """Get a query string for this reaction."""
        
        rdict = {-1: [], 1: []}
        for c_w_coeff in self.reactants:
            c = numpy.abs(c_w_coeff.coeff)
            s = numpy.sign(c_w_coeff.coeff)
            if s == 0:
                continue
            if c == 1:
                rdict[s].append(c_w_coeff.GetName())
            else:
                rdict[s].append('%d %s' % (c, c_w_coeff.GetName()))
                
        return '%s = %s' % (' + '.join(rdict[-1]), ' + '.join(rdict[1]))
    
    def ContainsCO2(self):
        return self._FindCompoundIndex('C00011') is not None
            
    def IsReactantFormulaMissing(self):
        for compound_w_coeff in self.reactants:
            if compound_w_coeff.compound.GetAtomBag() is None:
                return True
        
        return False

    def GetReactantFormulaMissing(self):
        for compound_w_coeff in self.reactants:
            if compound_w_coeff.compound.GetAtomBag() is None:
                yield compound_w_coeff
    
    def IsEmpty(self):
        return len(self.reactants) == 0
    
    def IsBalanced(self):
        """Checks if the collection is atom-wise balanced.
        
        Returns:
            True if the collection is atom-wise balanced.
        """
        try:
            return self._IsBalanced(self._GetAtomDiff())
        except ReactantFormulaMissingError:
            return False
    
    def IsElectronBalanced(self):
        """Checks if the collection is electron-wise balanced.
        
        Returns:
            True if the collection is electron-wise balanced.
        """
        return self._GetElectronDiff() == 0
    
    def StandardizeHalfReaction(self):
        """Checks if the reaction is a half-reaction (excess electrons).
        
        Returns:
            True if the collection is electron-wise balanced.
        """
        if self._GetElectronDiff() < 0:
            self.SwapSides()
            
    def E0_prime(self):
        """
            Returns the standard transformed reduction potential of this reaction.
        """
        delta_electrons = self._GetElectronDiff()
        assert delta_electrons != 0
        return - self.DeltaG0Prime() / (constants.F*delta_electrons)

    def Em_prime(self):
        """
            Returns the standard transformed reduction potential of this reaction.
        """
        delta_electrons = self._GetElectronDiff()
        assert delta_electrons != 0
        return - self.DeltaGmPrime() / (constants.F*delta_electrons)

    def E_prime(self):
        """Returns the standard transformed reduction potential of this reaction."""
        delta_electrons = self._GetElectronDiff()
        assert delta_electrons != 0
        return - self.DeltaGPrime() / (constants.F*delta_electrons)

    def EUncertainty(self):
        dg_u = self.DeltaGUncertainty()
        if dg_u is None:
            return None
        else:
            delta_electrons = self._GetElectronDiff()
            assert delta_electrons != 0
            return numpy.abs(dg_u / (constants.F*delta_electrons))
    
    def _ExtraWaters(self):
        atom_diff = self._GetAtomDiff()
            
        if not atom_diff:
            return None
                
        # Ignore hydrogen.
        atom_diff.pop('H', 0)
        
        # Omit oxygen for checking balancedness.
        oxy_count = atom_diff.pop('O', 0)
    
        # If it's not balanced without oxygen, can't balance with water.
        if not self._IsBalanced(atom_diff):
            return None
        
        # Requires this many waters to balance (1 O per).
        return oxy_count

    def _FindCompoundIndex(self, kegg_id):
        """
            Returns the first index of the compound with the given KEGG id.
        """
        logging.debug('Looking for the index of %s' % kegg_id)
        for i, c in enumerate(self.reactants):
            if c.compound.kegg_id == kegg_id:
                logging.debug('Found %s at index %d' % (kegg_id, i))
                return i
        return None

    def _AddCompound(self, kegg_id, how_many):
        """Adds "how_many" of the compound with the given id.
        
        Args:
            kegg_id: the KEGG id.
            how_many: by how much to change the reactant's coefficient.
        """
        i = self._FindCompoundIndex(kegg_id)
        if i is not None:
            self.reactants[i].coeff += how_many
        else:
            self.reactants += [CompoundWithCoeff.FromId(how_many, kegg_id)]
            
        # clear the cache since the reaction has changed
        self._catalyzing_enzymes = None

    def _ReplaceCompound(self, from_id, to_id):
        """
            Replace a compound in the reaction with another compound according
            to their IDs.
            The stoichiometric coefficient and concentration are copied from
            the old compound to the new one.
        """
        
        # set the coefficient of the original compound to 0
        i = self._FindCompoundIndex(from_id)
        if i is None:
            return
        how_many = self.reactants[i].coeff
        self.reactants[i].coeff = 0

        # create a new compound with the new kegg_id and the same coefficient
        # or add the number to the coefficient if it already is a reactant
        j = self._FindCompoundIndex(to_id)
        if j is None:
            self.reactants[i] = CompoundWithCoeff.FromId(how_many, to_id)
        else:
            self.reactants[j] += how_many
            self._Dedup()

        # clear the cache since the reaction has changed
        self._catalyzing_enzymes = None

    def _Dedup(self):
        """
            Collapses duplicate compounds and removes ones with coefficients of zero.
        """
        kegg_id_to_index = {}
        for i, c in enumerate(self.reactants):
            first_i = kegg_id_to_index.setdefault(c.compound.kegg_id, i)
            if i != first_i:
                self.reactants[first_i].coeff += c.coeff
                c.coeff = 0
                
        self.reactants = filter(lambda x: x.coeff != 0, self.reactants)
        
        # always make sure that H2O is the last reactant (so that it will
        # appear last in the chemical formula)
        i_h2o = self._FindCompoundIndex('C00001')
        if i_h2o is not None:
            self.reactants = self.reactants[:i_h2o] + \
                             self.reactants[(i_h2o + 1):] + \
                             [self.reactants[i_h2o]]

    def _FilterProtons(self):
        """
            Removes Protons from the list of compounds.
            Since we use Bob Alberty's framework for biochemical reactions, there
            is no meaning for having H+ in a reaction.
        """
        self.reactants = filter(lambda c: c.compound.kegg_id != 'C00080', self.reactants)
        
    def TryBalanceWithWater(self):
        """Try to balance the reaction with water.
        
        Returns:
            True if the reaction is balanced already or with
            additional waters on either side.
        """ 
        extra_waters = self._ExtraWaters()
        if extra_waters is None:
            # cannot balance the reaction with H2O only
            return False
        if extra_waters != 0:
            self._AddCompound('C00001', extra_waters)
            self._Dedup()
        return True
    
    def CanBalanceWithWater(self):
        """Returns True if balanced with or without water."""
        try:
            extra_waters = self._ExtraWaters()
        except ReactantFormulaMissingError:
            return True
        
        return extra_waters is not None
    
    def BalanceElectrons(self,
                         acceptor_id='C00003',           # NAD+
                         reduced_acceptor_id='C00004'):  # NADH
        """Try to balance the reaction electons."""        
        net_electrons = self._GetElectronDiff()
        if net_electrons != 0:
            self._AddCompound(reduced_acceptor_id, -net_electrons/2)
            self._AddCompound(acceptor_id, net_electrons/2)
            self._Dedup()
        
    def _GetConcentrationCorrection(self):
        """
            Get the concentration term in DeltaG' for these concentrations.
        
            Returns:
                The correction or None on error.
        """        
        # Shorthand for coeff * log(concentration)
        mult_log_c_list = [c.coeff * numpy.log(c.phase.Value())
                           for c in self.reactants]

        # Compute log(Q) - the log of the reaction quotient
        log_Q = sum(mult_log_c_list)
        
        _r = constants.R
        _t = constants.DEFAULT_TEMP
        return _r * _t * log_Q

    def _GetConcentrationCorrectionMilliMolar(self):
        """
            Get the concentration term in DeltaGm'.
                    
            Returns:
                The correction or None on error.
        """
        # calculate stoichiometric imbalance (how many more products are there
        # compared to substrates). Note that H2O isn't counted
        imbalance = sum([c.coeff for c in self.reactants
                         if not c.phase.IsConstant()])
        
        _r = constants.R
        _t = constants.DEFAULT_TEMP
        return _r * _t * numpy.log(1e-3) * imbalance

    def DeltaG0Prime(self):
        """Compute the DeltaG0' for a reaction.
        
        Returns:
            The DeltaG0' for this reaction, or None if data was missing.
        """
        if self._dg0_prime is not None:
            return self._dg0_prime

        logging.debug('Aqueous Params = ' + str(self.aq_params))
        c_dg0_prime_list = [c.DeltaG0Prime(self.aq_params) for c in self.reactants]

        # find all the IDs of compounds that have no known formation energy
        # if there are any such compounds, print and error message and 
        # return None since we cannot calculate the reaction energy
        kegg_id_list = [c.kegg_id for c in self.reactants]
        kegg_id_and_dg0 = zip(kegg_id_list, c_dg0_prime_list)
        unknown_kegg_ids = [x[0] for x in kegg_id_and_dg0 if x[1] is None]
        if unknown_kegg_ids:
            logging.warning("Failed to get formation energy for: " +
                            ', '.join(unknown_kegg_ids))
            return None
        self._dg0_prime = sum(c_dg0_prime_list)
        return self._dg0_prime

    def DeltaGmPrime(self):
        """Compute the DeltaGm' for a reaction (i.e. at 1 mM).
        
        Returns:
            The DeltaGm' for this reaction, or None if data was missing.
        """
        dg0_prime = self.DeltaG0Prime()
        correction = self._GetConcentrationCorrectionMilliMolar()
        return dg0_prime + correction

    def DeltaGPrime(self):
        """Compute the DeltaG' for a reaction.
        
        Returns:
            The DeltaG' for this reaction, or None if data was missing.
        """
        dg0_prime = self.DeltaG0Prime()
        correction = self._GetConcentrationCorrection()
        return dg0_prime + correction
        
    def HalfReactionDeltaGPrime(self):
        """Compute the DeltaG' for a half-reaction, assuming the missing
           electrons are provided in a certain potential 
           (e_reduction_potential)
        
        Returns:
            The DeltaG' for this half-reaction, or None if data was missing.
        """
        dg_prime = self.DeltaGPrime()
        delta_electrons = self._GetElectronDiff()  
        return dg_prime + constants.F * delta_electrons * \
                          self.aq_params.e_reduction_potential
    
    def KeqPrime(self):
        """
            Returns the transformed equilibrium constant for this reaction.
        """
        dg0_prime = self.DeltaG0Prime()
        if dg0_prime is None:
            return None
        
        rt = constants.R * constants.DEFAULT_TEMP
        keq = numpy.exp(-dg0_prime / rt)
        return keq

    def KeqPrimeHuman(self):
        """
            Returns the transformed equilibrium constant for this reaction,
            in a human readable formet (using HTML superscript).
        """
        dg0_prime = self.DeltaG0Prime()
        if dg0_prime is None:
            return None
        
        rtln10 = constants.R * constants.DEFAULT_TEMP * numpy.log(10)
        x = -dg0_prime / rtln10

        exp = numpy.floor(x)
        prefactor = 10**(x - exp)
        if abs(exp) <= 2:
            return '%.2g' % (10**x)
        else:
            return '%.1f &times; 10<sup>%d</sup>' % (prefactor, exp)
        
    def NoDeltaGExplanation(self):
        """
            Get an explanation for why there's no delta G value.
        
            Return:
                The explanation or None.
        """
        for compound in self.reactants:
            if compound.compound.no_dg_explanation:
                name = compound.compound.common_names.all()[0].name
                return '%s %s' % (name,
                                  compound.compound.no_dg_explanation.lower())
        return None
        
    def DeltaGUncertainty(self):
        if self._GetMaxCommonPriority() != 1:
            return None
        C1 = cc_preprocess['C1']
        C2 = cc_preprocess['C2']
        C3 = cc_preprocess['C3']

        Nc = C1.shape[0]        
        Ng = C3.shape[0]      

        assert C1.shape[0] == C1.shape[1]
        assert C1.shape[1] == C2.shape[0]
        assert C2.shape[1] == C3.shape[0]
        assert C3.shape[0] == C3.shape[1]
        
        # x is the stoichiometric vector of the reaction, only for the
        # compounds that appeared in the original training set for CC
        x = numpy.matrix(numpy.zeros((Nc, 1))) 

        # g is the group incidence vector of all the other compounds
        g = numpy.matrix(numpy.zeros((Ng, 1)))
        logging.debug('g.shape = %s' % str(g.shape))
        for compound in self.reactants:
            logging.debug(compound.compound.kegg_id)
            i = compound.compound.index
            gv = compound.compound.sparse_gv
            if i is not None:
                logging.debug('index = %d' % i)
                x[i, 0] = compound.coeff
            elif gv is not None:
                logging.debug('len(gv) = %d' % len(gv))
                for g_ind, g_count in gv:
                    g[g_ind, 0] += g_count
            else:
                logging.debug('could not find index nor group vector')
                return None

        s_cc = float(numpy.sqrt(x.T * C1 * x + x.T * C2 * g + g.T * C3 * g))
        logging.debug('s_cc = %g' % s_cc)
        return 1.96*s_cc

    def AllCompoundsWithTransformedEnergies(self):
        for c_w_coeff in self.reactants:
            dgt = c_w_coeff.compound.DeltaG0Prime()
            c_w_coeff.transformed_energy = dgt
            yield c_w_coeff

    def ExtraAtoms(self):
        try:
            diff = self._GetAtomDiff()
        except ReactantFormulaMissingError:
            return None
        diff.pop('H', 0)
        extras = filter(lambda t: t[1] > 0, diff.iteritems())
        if not extras:
            return None
        
        extras.sort(key=lambda t: t[1], reverse=True)
        return extras

    def MissingAtoms(self):
        try:
            diff = self._GetAtomDiff()
        except ReactantFormulaMissingError:
            return None
        diff.pop('H', 0)
        short = filter(lambda t: t[1] < 0, diff.iteritems())
        if not short:
            return None
        
        short = [(atom, -count) for atom, count in short]
        short.sort(key=lambda t: t[1], reverse=True)        
        return short
    
    def ExtraElectrons(self):
        diff = self._GetElectronDiff()
        if diff > 0:
            return diff
        return None
    
    def MissingElectrons(self):
        diff = self._GetElectronDiff()
        if diff < 0:
            return -diff
        return None
    
    def _CheckConservationLaw(self, sparse_reaction, claw):
        inner_prod = 0.0
        for kegg_id, coeff in claw.GetSparseRepresentation().iteritems():
            if kegg_id in sparse_reaction:
                inner_prod += coeff * sparse_reaction[kegg_id]
        return abs(inner_prod) < 1e-10
    
    def CheckConservationLaws(self):
        sparse_reaction = dict([(c.compound.kegg_id, c.coeff) 
                                for c in self.reactants])

        all_claws = models.ConservationLaw.objects.select_related().all()
        for claw in all_claws:
            if not self._CheckConservationLaw(sparse_reaction, claw):
                return False
        return True
        
    def _IsPhysiologicalConcentration(self):
        return all(c_w_c.phase.IsPhysiological() for c_w_c in self.reactants)
    
    def GetSourceReference(self):
        """
            Assuming that all reactants in the chosen priority group have
            the same source reference, we provide the ref from the first
            reactant to represent all of them.
        """
        #return self.reactants[0].compound._GetDGSource()
        source_names = set(map(lambda x : str(x.compound._GetDGSource()), self.reactants))
        return ', '.join([self.GetSourceReferenceLink(s) for s in source_names])
    
    def GetSourceReferenceLink(self, source_name):
        try:
            source = models.ValueSource.objects.get(name=source_name)
            url = source.url
        except Exception:
            url = None
        
        if url:
            return '<a href="%s">%s</a>' % (url, source_name)
        else:
            return '<a href="/data_refs">%s</a>' % (source_name)
        
    
    substrates = property(GetSubstrates)
    products = property(GetProducts)
    contains_co2 = property(ContainsCO2)
    is_conserving = property(CheckConservationLaws)
    is_reactant_formula_missing = property(IsReactantFormulaMissing)
    reactants_with_missing_formula = property(GetReactantFormulaMissing)
    is_empty = property(IsEmpty)    
    is_balanced = property(IsBalanced)
    is_electron_balanced = property(IsElectronBalanced)
    balanced_with_water = property(CanBalanceWithWater)
    extra_atoms = property(ExtraAtoms)
    missing_atoms = property(MissingAtoms)
    extra_electrons = property(ExtraElectrons)
    missing_electrons = property(MissingElectrons)
    all_compounds = property(AllCompoundsWithTransformedEnergies)
    dg0_prime = property(DeltaG0Prime)
    dgm_prime = property(DeltaGmPrime)
    dg_prime = property(DeltaGPrime)
    half_reaction_dg_prime = property(HalfReactionDeltaGPrime)
    k_eq_prime = property(KeqPrime)
    k_eq_prime_human = property(KeqPrimeHuman)
    e0_prime = property(E0_prime)
    em_prime = property(Em_prime)
    e_prime = property(E_prime)
    e_uncertainty = property(EUncertainty)
    no_dg_explanation = property(NoDeltaGExplanation)
    dg_uncertainty = property(DeltaGUncertainty)
    cc_link = property(lambda self: self.GetNewPriorityLink(1))
    alberty_link = property(lambda self: self.GetNewPriorityLink(99))
    ph_graph_link = property(GetPhGraphLink)
    pmg_graph_link = property(GetPMgGraphLink)
    is_graph_link = property(GetIonicStrengthGraphLink)
    catalyzing_enzymes = property(_GetCatalyzingEnzymes)
    is_phys_conc = property(_IsPhysiologicalConcentration)
    source_reference = property(GetSourceReference)
