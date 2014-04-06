# -*- coding: utf-8 -*-

import logging
import numpy
import urllib

from gibbs import conditions
from gibbs import constants
from gibbs import models

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
        self._name = name
        self.transformed_energy = None
        self.phase = phase
    
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
    def FromId(coeff, kegg_id, name=None):
        try:
            compound = models.Compound.objects.get(kegg_id=kegg_id)
        except Exception:
            return None
        my_name = name or compound.FirstName()
        
        return CompoundWithCoeff(coeff, compound,
                                 phase=None, name=my_name)
        
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

    def DeltaG0Prime(self, pH, pMg, ionic_strength):
        dg0_prime = self.compound.DeltaG0Prime(pH=pH,
                                               pMg=pMg,
                                               ionic_strength=ionic_strength,
                                               phase=self.phase.PhaseName())
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
        return self.phase.HumanValueAndUnits()[0]

    def GetPhasePrefactor(self):
        return self.phase.HumanValueAndUnits()[1]
        
    def GetPhaseUnits(self):
        return self.phase.HumanValueAndUnits()[2]

    def GetPhaseIsConstant(self):
        return self.phase.IsConstant()
        
    def GetPhaseHumanString(self):
        return str(self.phase)
    
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
    
    def __init__(self, reactants=None,
                 pH=constants.DEFAULT_PH,
                 pMg=constants.DEFAULT_PMG,
                 ionic_strength=constants.DEFAULT_IONIC_STRENGTH,
                 e_reduction_potential=constants.DEFAULT_ELECTRON_REDUCTION_POTENTIAL):
        """Construction.
        
        Args:
            substrates: a list of CompoundWithCoeff objects.
            products: a list of CompoundWithCoeff objects.
        """
        self.reactants = reactants or []
        self._FilterProtons()
        self._Dedup()
        
        self.ph = pH
        self.pmg = pMg
        self.i_s = ionic_strength
        self.e_reduction_potential = e_reduction_potential
        self._dg0_prime = None
        self._conditions = None
        self._kegg_id = None

        # used only as cache, no need to copy while cloning        
        self._catalyzing_enzymes = None
        
        self._SetCompoundPriorities()

    def Clone(self):
        """
            Make a copy of the 
        """
        rs = [r.Clone() for r in self.reactants]
        other = Reaction(rs, self.ph, self.pmg, self.i_s,
                         self.e_reduction_potential)
        other._dg0_prime = self._dg0_prime
        other._conditions = self._conditions
        return other

    def _SetCompoundPriorities(self):
        """
            Returns a set of (int, SpeciesGroup) tuples for the reaction.
        """
        
        priorities = [c.compound.GetSpeciesGroupPriorities() for c in self.reactants]
        priorities = filter(lambda l: len(l) > 0, priorities)
        
        # Someone is missing data!
        if priorities == []:
            return
        priority_to_use = min([max(l) for l in priorities if l])
        logging.debug('All priorities: %s' % str(priorities))        
        logging.debug('Using the priority %d' % priority_to_use)        
        
        for c in self.reactants:
            c.compound.SetSpeciesGroupPriority(priority_to_use)
         
    def SwapSides(self):
        """Swap the sides of this reaction."""
        for c in self.reactants:
            c.coeff = -c.coeff
        
    def GetConditions(self):
        """Get the concentration profile of this reaction."""
        return self._conditions
    
    def GetSubstrates(self):
        s = [c for c in self.reactants if c.coeff < 0]
        return s
        
    def GetProducts(self):
        p = [c for c in self.reactants if c.coeff > 0]
        return p
    
    def SetConditions(self, cond):
        """Apply this concentration profile to this reaction.
        
        Args:
            cond: a _BaseConditions object.
        """
        self._dg0_prime = None
        self._conditions = cond

        for c in self.reactants:
            c.phase = self._conditions.GetPhase(c.kegg_id)
            if c.phase is None:
                phase_name = c.compound.GetDefaultPhaseName()
                self._conditions.SetPhase(c.kegg_id, phase_name)
                c.phase = self._conditions.GetPhase(c.kegg_id)
                logging.info('Using default phase for %s: %s' %
                             (c.kegg_id, phase_name))

    conditions = property(GetConditions, SetConditions)
    
    def __str__(self):
        """
            Simple text reaction representation.
        """
        rlist = map(str, self.substrates)
        plist = map(str, self.products)
        return '%s <=> %s' % (' + '.join(rlist), ' + '.join(plist))

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
                'pH': self.ph,
                'ionic_strength': self.i_s}
        if self.k_eq_prime is not None:
            d['keq_prime'] = {
                'value': self.k_eq_prime,
                'pH': self.ph,
                'ionic_strength': self.i_s}
        
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
        if form.cleaned_reactionId:
            stored_reaction = models.StoredReaction.objects.get(
                kegg_id=form.cleaned_reactionId)
            return Reaction.FromStoredReaction(stored_reaction)

        ph = form.cleaned_ph
        pmg = form.cleaned_pmg
        i_s = form.cleaned_ionic_strength
        e_red = form.cleaned_e_reduction_potential
        logging.info('pH = %.1f, pMg = %.1f, ionic strength = %.1f, Ered = %.1f' %
                     (ph, pmg, i_s, e_red))
        
        if form.cleaned_reactantsPhase == []:
            phases = [None] * len(form.cleaned_reactantsCoeff)
        else:
            phases = form.cleaned_reactantsPhase

        zipped_reactant_data = zip(form.cleaned_reactantsCoeff,
                                   form.cleaned_reactantsId,
                                   phases,
                                   form.cleaned_reactantsName)

        compound_list = []
        for coeff, kegg_id, phase, name in zipped_reactant_data:
            compound_list.append({'coeff': coeff,
                                  'kegg_id': kegg_id,
                                  'phase': phase,
                                  'name': name})

        # Build the appropriate concentration profile.
        cond = conditions.CreateConditions(form.cleaned_conditions,
                                           form.cleaned_reactantsId,
                                           form.cleaned_reactantsPhase,
                                           form.cleaned_reactantsConcentration)
        
        # Return the built reaction object.
        return Reaction.FromIds(compound_list,
                                cond=cond,
                                pH=ph, pMg=pmg,
                                ionic_strength=i_s,
                                e_reduction_potential=e_red)
    
    @staticmethod
    def FromIds(compound_list,
                cond=None,
                pH=constants.DEFAULT_PH,
                pMg=constants.DEFAULT_PMG,
                ionic_strength=constants.DEFAULT_IONIC_STRENGTH,
                e_reduction_potential=constants.DEFAULT_ELECTRON_REDUCTION_POTENTIAL):
        """Build a reaction object from lists of IDs.
        
        Args:
            compound_list: an iterable of (coeff, kegg_id, name) of reactants.
            concentration_profile: a ConcentrationProfile object.
            
        Returns:
            A properly set-up Reaction object or None if there's an error.
        """
        cond = cond or conditions.StandardConditions()
        
        reactants = []
        # Get products and substrates.
        for c in compound_list:
            coeff = c['coeff']
            kegg_id = c['kegg_id']
            phase = c.get('phase', None)
            name = c.get('name', None)
            logging.info('Adding compound %s with coeff %d and phase %s' %
                         (kegg_id, coeff, phase))
            c_w_c = CompoundWithCoeff.FromId(coeff, kegg_id, name=name)
            reactants.append(c_w_c)
            
            if phase is not None:
                cond.SetPhase(kegg_id, phase)

        rxn = Reaction(reactants, pH=pH, pMg=pMg,
                       ionic_strength=ionic_strength,
                       e_reduction_potential=e_reduction_potential)
        
        rxn.conditions = cond or conditions.StandardConditions()
        
        return rxn
    
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
        """Get the URL params for this reaction."""
        params = []
        for compound in self.reactants:
            kegg_id = compound.compound.kegg_id
            params.append('reactantsId=%s' % kegg_id)
            params.append('reactantsCoeff=%d' % compound.coeff)
            if compound.name:
                params.append('reactantsName=%s' % compound.name)
            if self._conditions:
                params.extend(self._conditions._GetUrlParams(kegg_id))
        
        if self.ph:
            params.append('ph=%f' % self.ph)
        if self.pmg:
            params.append('pmg=%f' % self.pmg)
        if self.i_s:
            params.append('ionic_strength=%f' % self.i_s)
                
        if query:
            for arrow in constants.POSSIBLE_REACTION_ARROWS:
                tmp_query = query.replace(arrow, '=>')
            params.append('query=%s' % urllib.quote(tmp_query))
            
        return params
    
    def GetHyperlink(self, query=None):
        return '/reaction?%s' % '&'.join(self._GetUrlParams(query))
    
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
        return other.GetHyperlink(query)

    def GetHalfReactionLink(self, query=None):
        """
            Returns a link to the same reaction,
            balanced by extra H2O molecules.
        """
        params = self._GetUrlParams(query)
        return '/half_reaction?%s' % '&'.join(params)

    def GetTemplateData(self, query=None):
        template_data = {'reaction': self,
                         'query': query}
        try:
            balance_with_water_link = self.GetBalanceWithWaterLink(query)
            balance_electrons_link = self.GetBalanceElectronsLink(query)
            half_reaction_link = self.GetHalfReactionLink(query)
            replace_co2_link = self.GetReplaceCO2Link(query)
            template_data.update({
                'balance_with_water_link': balance_with_water_link,
                'balance_electrons_link': balance_electrons_link,
                'replace_co2_link': replace_co2_link,
                'half_reaction_link': half_reaction_link})
        except ReactantFormulaMissingError:
            pass
        if self._conditions is not None:
            template_data.update(self._conditions.GetTemplateDict())
        return template_data
    
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
                
        return '%s <=> %s' % (' + '.join(rdict[-1]), ' + '.join(rdict[1]))
    
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
            
    def E0_prime(self, pH=None, pMg=None, ionic_strength=None):
        """Returns the standard transformed reduction potential of this reaction."""
        delta_electrons = self._GetElectronDiff()
        assert delta_electrons != 0
        return - self.DeltaG0Prime(pH, pMg, ionic_strength) / (constants.F*delta_electrons)

    def E_prime(self, pH=None, pMg=None, ionic_strength=None):
        """Returns the standard transformed reduction potential of this reaction."""
        delta_electrons = self._GetElectronDiff()
        assert delta_electrons != 0
        return - self.DeltaGPrime(pH, pMg, ionic_strength) / (constants.F*delta_electrons)
    
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
            c_w_coeff = CompoundWithCoeff.FromId(how_many, kegg_id)
            self.reactants.append(c_w_coeff)
            self._conditions.SetPhase(kegg_id,
                                      c_w_coeff.compound.GetDefaultPhaseName())

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
            c_w_c = CompoundWithCoeff.FromId(how_many, to_id)
            self._conditions.SetPhase(to_id,
                                      c_w_c.compound.GetDefaultPhaseName())
            self.reactants[i] = c_w_c
        else:
            self.reactants[j] += how_many
            self._Dedup()

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
        """Get the concentration term in DeltaG' for these concentrations.
        
        Args:
            concentration_profile: a ConcentrationProfile object.
        
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

    def DeltaG0Prime(self, pH=None, pMg=None, ionic_strength=None):
        """Compute the DeltaG0' for a reaction.
        
        Returns:
            The DeltaG0' for this reaction, or None if data was missing.
        """
        ph = pH or self.ph
        pmg = pMg or self.pmg
        i_s = ionic_strength or self.i_s

        if self._dg0_prime is not None:
            if (ph, pmg, i_s) == (self.ph, self.pmg, self.i_s):
                return self._dg0_prime
            else:
                self._dg0_prime = None

        c_dg0_prime_list = [c.DeltaG0Prime(pH=ph, pMg=pmg, ionic_strength=i_s)
                            for c in self.reactants]

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

    def DeltaGPrime(self, pH=None, pMg=None, ionic_strength=None):
        """Compute the DeltaG' for a reaction.
        
        Returns:
            The DeltaG' for this reaction, or None if data was missing.
        """
        dg0_prime = self.DeltaG0Prime(pH=pH, pMg=pMg,
                                      ionic_strength=ionic_strength)
        correction = self._GetConcentrationCorrection()
        return dg0_prime + correction
        
    def HalfReactionDeltaGPrime(self, pH=None, pMg=None, ionic_strength=None,
                                e_reduction_potential=None):
        """Compute the DeltaG' for a half-reaction, assuming the missing
           electrons are provided in a certain potential 
           (e_reduction_potential)
        
        Returns:
            The DeltaG' for this half-reaction, or None if data was missing.
        """
        dg_prime = self.DeltaGPrime(pH=pH, pMg=pMg,
                                    ionic_strength=ionic_strength)
        e_red = e_reduction_potential or self.e_reduction_potential
        delta_electrons = abs(self._GetElectronDiff())      
        return dg_prime + constants.F * delta_electrons * e_red
    
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

    def AllCompoundsWithTransformedEnergies(self):
        for c_w_coeff in self.reactants:
            dgt = c_w_coeff.compound.DeltaG0Prime(pH=self.ph,
                                                  pMg=self.pmg,
                                                  ionic_strength=self.i_s)
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
    dg_prime = property(DeltaGPrime)
    half_reaction_dg_prime = property(HalfReactionDeltaGPrime)
    k_eq_prime = property(KeqPrime)
    k_eq_prime_human = property(KeqPrimeHuman)
    e0_prime = property(E0_prime)
    e_prime = property(E_prime)
    no_dg_explanation = property(NoDeltaGExplanation)
    ph_graph_link = property(GetPhGraphLink)
    pmg_graph_link = property(GetPMgGraphLink)
    is_graph_link = property(GetIonicStrengthGraphLink)
    catalyzing_enzymes = property(_GetCatalyzingEnzymes)
