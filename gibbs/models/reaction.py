# -*- coding: utf-8 -*-

import hashlib
import logging
import numpy
import urllib
import os
import json
from scipy.sparse import csr_matrix
from django.db import models
from django.apps import apps
from util import constants
from .. import conditions
from .compound import CommonName, CompoundWithCoeff


class Preprocessing(object):
    relpath = os.path.dirname(os.path.realpath(__file__))
    cc_preprocess_fname = os.path.join(relpath, '../../data/cc_preprocess.npz')
    cc_preprocess = numpy.load(cc_preprocess_fname)

    C1 = numpy.matrix(cc_preprocess['C1'])
    C2 = numpy.matrix(cc_preprocess['C2'])
    C3 = numpy.matrix(cc_preprocess['C3'])
    G1 = numpy.matrix(cc_preprocess['G1'])
    G2 = numpy.matrix(cc_preprocess['G2'])
    G3 = numpy.matrix(cc_preprocess['G3'])
    S = numpy.matrix(cc_preprocess['S'])
    cids = cc_preprocess['cids']
    Nc = C1.shape[0]
    Ng = C3.shape[0]
    assert C1.shape[0] == C1.shape[1]
    assert C1.shape[1] == C2.shape[0]
    assert C2.shape[1] == C3.shape[0]
    assert C3.shape[0] == C3.shape[1]
    assert C3.shape[0] == C3.shape[1]

    @staticmethod
    def GetCompoundVectors(compound):
        # x is the stoichiometric vector of the reaction, only for the
        # compounds that appeared in the original training set for CC
        x = numpy.matrix(numpy.zeros((Preprocessing.Nc, 1)))

        # g is the group incidence vector of all the other compounds
        g = numpy.matrix(numpy.zeros((Preprocessing.Ng, 1)))
        logging.debug(compound.compound.kegg_id)
        i = compound.compound.index
        gv = compound.compound.sparse_gv
        if i is not None:
            logging.debug('index = %d' % i)
            x[i, 0] = compound.coeff
        elif gv is not None:
            logging.debug('len(gv) = %d' % len(gv))
            for g_ind, g_count in gv:
                g[g_ind, 0] += compound.coeff * g_count
        else:
            raise Exception('could not find index nor group vector for %s'
                            % compound.compound.kegg_id)
        return x, g

    @staticmethod
    def GetReactionVectors(reactants):
        # x is the stoichiometric vector of the reaction, only for the
        # compounds that appeared in the original training set for CC
        x_reaction = numpy.matrix(numpy.zeros((Preprocessing.Nc, 1)))

        # g is the group incidence vector of all the other compounds
        g_reaction = numpy.matrix(numpy.zeros((Preprocessing.Ng, 1)))
        for x, g in map(Preprocessing.GetCompoundVectors, reactants):
            x_reaction += x
            g_reaction += g

        logging.debug('x = %s' % csr_matrix(x_reaction))
        logging.debug('g = %s' % csr_matrix(g_reaction))
        return x_reaction, g_reaction

    @staticmethod
    def DeltaGUncertainty(x, g):
        return float(numpy.sqrt(x.T * Preprocessing.C1 * x +
                             x.T * Preprocessing.C2 * g +
                             g.T * Preprocessing.C3 * g))

    @staticmethod
    def WriteCompoundAndCoeff(kegg_id, coeff):
        if coeff == 1:
            return kegg_id
        else:
            return "%g %s" % (coeff, kegg_id)

    @staticmethod
    def DictToReactionString(d):
        """String representation."""
        left = []
        right = []
        for kegg_id, coeff in sorted(d.items()):
            _s = Preprocessing.WriteCompoundAndCoeff(kegg_id, -coeff)
            if coeff < 0:
                left.append(_s)
            elif coeff > 0:
                right.append(_s)
        return "%s %s %s" % (' + '.join(left), '<=>', ' + '.join(right))

    @staticmethod
    def Analyze(x, g):
        weights_rc = x.T * Preprocessing.G1
        weights_gc = x.T * Preprocessing.G2 + g.T * Preprocessing.G3
        weights = weights_rc + weights_gc

        res = []
        for j in range(Preprocessing.S.shape[1]):
            d = {Preprocessing.cids[i]: Preprocessing.S[i, j]
                 for i in range(Preprocessing.Nc)
                 if Preprocessing.S[i, j] != 0}
            r_string = Preprocessing.DictToReactionString(d)
            res.append({'w': weights[0, j],
                        'w_rc': weights_rc[0, j].round(4),
                        'w_gc': weights_gc[0, j].round(4),
                        'reaction_string': r_string})
        res.sort(key=lambda d: abs(d['w']), reverse=True)
        return res

    @staticmethod
    def IsUsingGroupContributions(x, g):
        weights_gc = x.T * Preprocessing.G2 + g.T * Preprocessing.G3
        sum_w_gc = sum(numpy.abs(weights_gc).flat)
        logging.debug('sum(w_gc) = %.2g' % sum_w_gc)
        return sum_w_gc > 1e-5


class ReactantFormulaMissingError(Exception):

    def __init__(self, c):
        self.compound = c

    def __str__(self):
        return ("Cannot test reaction balancing because the reactant"
                " %s does not have a chemical formula" % self.compound.kegg_id)


class Reaction(models.Model):
    """A reaction."""

    def __init__(self, reactants=None, aq_params=None):
        """Construction.

        Args:
            reactants: a list of CompoundWithCoeff objects.
        """
        self.reactants = reactants or []
        self.aq_params = aq_params or conditions.AqueousParams()

        self._is_formation_reaction = len(self.reactants) < 2

        if not self._is_formation_reaction:
            self._FilterProtonsAndElectrons()
            self._Dedup()

        self._kegg_id = None
        self._uncertainty = None
        self.is_using_gc = False

        # used only as cache, no need to copy while cloning
        self._catalyzing_enzymes = None

    def SetAqueousParams(self, aq_params):
        self._dg0_prime = None
        self._aq_params = aq_params
        self._SetCompoundPriorities()

    def GetAqueousParams(self):
        return self._aq_params

    aq_params = property(GetAqueousParams, SetAqueousParams)

    def ResetConcentrations(self):
        for c in self.reactants:
            c.ResetConcentration()

    def Clone(self):
        """
            Make a copy of the reaction
        """
        logging.debug('Cloning reaction...')
        other = Reaction()
        other.reactants = [r.Clone() for r in self.reactants]
        other._dg0_prime = self._dg0_prime
        other._uncertainty = self._uncertainty
        other.is_using_gc = self.is_using_gc
        other._aq_params = self._aq_params.Clone()
        other._kegg_id = self._kegg_id
        return other

    def _GetMaxCommonPriority(self):
        max_priority = self.aq_params.max_priority or 1

        # The chosen priority will be the highest number which is common
        # to all reactants.
        priorities = [c.compound.GetSpeciesGroupPriorities()
                      for c in self.reactants]
        priorities = list(filter(lambda l: len(l) > 0, priorities))

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

    @staticmethod
    def _GetStoredHashString(s):
        return apps.get_model('gibbs.StoredReaction').HashableReactionString(s)

    def GetHashableReactionString(self):
        return Reaction._GetStoredHashString(self.GetSparseRepresentation())

    def GetHash(self):
        return Reaction._GetStoredHashString(self.GetSparseRepresentation())

    def GetSpecialReactionWarning(self):

        def GetLearnMoreLink(faq_mark):
            return '</br><a href="/static/classic_rxns/faq.html#%s">Learn more &raquo;</a>' % faq_mark

        my_hash = self.GetHashableReactionString()

        atp_sparse = {'C00002': -1, 'C00001': -1, 'C00008': 1, 'C00009': 1}
        co2_sparse = {'C00011': -1, 'C00001': -1, 'C01353': 1}
        atp_hash = Reaction._GetStoredHashString(atp_sparse)
        co2_hash = Reaction._GetStoredHashString(co2_sparse)

        if my_hash == atp_hash:
            return ("The &Delta;G' of ATP hydrolysis is highly affected " +
                    "by Mg<sup>2+</sup> ions." +
                    GetLearnMoreLink("atp-hydrolysis"))
        elif my_hash == co2_hash:
            return ("You are looking at the &Delta;G' of CO<sub>2</sub> hydration." +
                    GetLearnMoreLink("co2-total"))
        elif (self._FindCompoundIndex('C14818') is not None and
              self._FindCompoundIndex('C14819') is not None):
            return ("Energetics of iron redox reactions depend heavily on the " +
                    "chemical forms of iron involved." +
                    GetLearnMoreLink("iron-redox"))
        elif (self._FindCompoundIndex('C00011') is not None and
              self._FindCompoundIndex('C00288') is None):
            return ('Did you mean <a href="%s">CO<sub>2</sub>(total)</a>?' % self.GetReplaceCO2Link() +
                    GetLearnMoreLink("co2-total"))
        elif (self._FindCompoundIndex('C00011') is not None and
              self._FindCompoundIndex('C00288') is not None):
            return ("One should not use CO<sub>2</sub>(aq) together with " +
                    "CO<sub>2</sub>(total) in the same reaction." +
                    GetLearnMoreLink("co2-total"))
        elif (self._FindCompoundIndex('C01353') is not None and
              self._FindCompoundIndex('C00288') is not None):
            return ("One should not use HCO<sub>3</sub><sup>-</sup>(aq) together with " +
                    "CO<sub>2</sub>(total) in the same reaction." +
                    GetLearnMoreLink("co2-total"))
        else:
            return False

    def _GetAllStoredReactions(self):
        """
            Find all stored reactions matching this compound (using the hash)
        """
        logging.debug('looking for stored reactions matching this one')
        my_hash = self.GetHash()
        my_string = self.GetHashableReactionString()

        matching_stored_reactions = apps.get_model(
            'gibbs.StoredReaction').objects.select_related().filter(
            reaction_hash=my_hash)
        logging.debug('my hash = %s (%d matches)' %
                      (my_hash, len(matching_stored_reactions)))

        matching_stored_reactions = \
            [m for m in matching_stored_reactions if
             m.GetHashableReactionString() == my_string]
        logging.debug('my hashable string = %s (%d matches)' %
                      (my_string, len(matching_stored_reactions)))

        return matching_stored_reactions

    @property
    def stored_reaction_id(self):
        stored_rxns = self._GetAllStoredReactions()
        if not stored_rxns:
            return None
        return stored_rxns[0].kegg_id

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
    def FromForm(form, aq_params=None):
        """
            Build a reaction object from a ReactionForm.

            Args:
                form: a ReactionForm object.

            Returns:
                A Reaction object or None if there's an error.
        """
#        max_priority = form.cleaned_max_priority
        n_react = len(form.cleaned_reactantsCoeff)
        coeffs = list(form.cleaned_reactantsCoeff)
        kegg_ids = list(form.cleaned_reactantsId)
        names = list(form.cleaned_reactantsName)
        phases = list(form.cleaned_reactantsPhase)
        concentrations = list(form.cleaned_reactantsConcentration)

        compound_list = []
        for i in range(n_react):
            d = {'coeff': coeffs[i], 'kegg_id': kegg_ids[i], 'name': names[i]}
            if phases != []:
                d['phase'] = phases[i]
            if concentrations != []:
                logging.debug("concentrations = %s" % str(concentrations))
                d['conc'] = concentrations[i]
            compound_list.append(d)

        # Return the built reaction object.
        return Reaction.FromIds(compound_list, aq_params=aq_params)

    @staticmethod
    def FromIds(compound_list, aq_params=None, fetch_db_names=False):
        """Build a reaction object from lists of IDs.

        Args:
            compound_list:  an iterable of dictionaries of reactants.
            keys:           kegg_id, coeff, phase, name
            aq_params:      an aqueous params object.
            fetch_db_names: if compound names should be fetched from the
                            database.

        Returns:
            A properly set-up Reaction object or None if there's an error.
        """
        kegg_ids = [d['kegg_id'] for d in compound_list]
        comps = apps.get_model('gibbs.Compound').objects.prefetch_related(
            'species_groups', 'species_groups__species',
            'common_names').filter(kegg_id__in=kegg_ids)
        kegg_id_to_compound = {c.kegg_id: c for c in comps}
        for d in compound_list:
            d['compound'] = kegg_id_to_compound[d['kegg_id']]
        if fetch_db_names:
            for d in compound_list:
                d['name'] = d['compound'].FirstName()
        reactants = list(map(CompoundWithCoeff.FromDict, compound_list))
        return Reaction(reactants, aq_params=aq_params)

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

            for atomic_number, atom_count in atom_bag.items():
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
            if electrons is None:
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
            params.append('query=%s' % urllib.parse.quote(tmp_query))

        return params

    def GetHyperlink(self, query=None):
        params = self._GetUrlParams(query)
        return '/reaction?%s' % '&'.join(params)

    def GetBalanceWithWaterLink(self, query=None):
        """Returns a link to balance this reaction with water."""
        other = self.Clone()
        other.TryBalanceWithWater()
        return other.GetHyperlink(query)

    def GetBalanceWithCoALink(self, query=None):
        """Returns a link to balance this reaction with water."""
        other = self.Clone()
        other.TryBalanceWithCoA()
        return other.GetHyperlink(query)

    def GetBalanceWithPiLink(self, query=None):
        """Returns a link to balance this reaction with water."""
        other = self.Clone()
        other.TryBalanceWithPi()
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
                         'query': query,
                         'alberty_link': self.GetNewPriorityLink(99),
                         'cc_link': self.GetNewPriorityLink(1),
                         'balance_with_water_link': None,
                         'balance_electrons_link': None,
                         'balance_with_coa_link': None,
                         'balance_with_pi_link': None}
        if not self._is_formation_reaction:
            try:
                if not self.IsBalanced():
                    template_data.update({'balance_with_water_link':
                                          self.GetBalanceWithWaterLink(query),
                                          'balance_with_coa_link':
                                          self.GetBalanceWithCoALink(query),
                                          'balance_with_pi_link':
                                          self.GetBalanceWithPiLink(query)})
                if not self.IsElectronBalanced():
                    template_data.update({'balance_electrons_link':
                                         self.GetBalanceElectronsLink(query)})
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
                rdict[s].append('%g %s' % (c, c_w_coeff.GetName()))

        return '%s <=> %s' % (' + '.join(rdict[-1]), ' + '.join(rdict[1]))

    def GetSlugQueryString(self):
        """Get an easily parsed query string for this reaction.

        This string contains no spaces or dashes in the compound names.
        Much easier to parse with REGEX. Used for generating SBtab.

        TODO unify with above?
        """
        rdict = {-1: [], 1: []}
        for c_w_coeff in self.reactants:
            c = numpy.abs(c_w_coeff.coeff)
            s = numpy.sign(c_w_coeff.coeff)
            if s == 0:
                continue
            if c == 1:
                rdict[s].append(c_w_coeff.compound.name_slug)
            else:
                rdict[s].append('%g %s' % (c, c_w_coeff.compound.name_slug))

        return '%s <=> %s' % (' + '.join(rdict[-1]), ' + '.join(rdict[1]))

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
            Returns the standard transformed reduction potential
            of this reaction.
        """
        delta_electrons = self._GetElectronDiff()
        assert delta_electrons != 0
        return - self.DeltaG0Prime() / (constants.F*delta_electrons)

    def Em_prime(self):
        """
            Returns the standard transformed reduction potential
            of this reaction.
        """
        delta_electrons = self._GetElectronDiff()
        assert delta_electrons != 0
        return - self.DeltaGmPrime() / (constants.F*delta_electrons)

    def E_prime(self):
        """
            Returns the standard transformed reduction
            potential of this reaction.
        """
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

    def _ExtraPis(self):
        # TODO probably could write a generic version of this.
        atom_diff = self._GetAtomDiff()

        if not atom_diff:
            return None

        # Pi = P O4
        # Ignore hydrogen as usual.
        atom_diff.pop('H', 0)

        os = atom_diff.pop('O', 0)
        ps = atom_diff.pop('P', 0)

        n_pis = numpy.array([ps / 1.0, os / 4.0])

        # Don't need any other elements
        pi_completes = self._IsBalanced(atom_diff)
        # Stoichiometry right for pi
        pi_completes &= numpy.all(n_pis == n_pis.astype(numpy.int))
        pi_completes &= numpy.all(n_pis ==
                               numpy.ones(n_pis.size) * n_pis[0])

        extra_pis = int(n_pis[0])
        if not pi_completes or not extra_pis:
            return None

        # Number of missing Pis
        return extra_pis

    def _ExtraCoAs(self):
        atom_diff = self._GetAtomDiff()

        if not atom_diff:
            return None

        # CoA = C21 H36 N7 016 P3 S
        # Ignore hydrogen as usual.
        atom_diff.pop('H', 0)

        cs = atom_diff.pop('C', 0)
        os = atom_diff.pop('O', 0)
        ns = atom_diff.pop('N', 0)
        ps = atom_diff.pop('P', 0)
        ss = atom_diff.pop('S', 0)

        n_coas = numpy.array([cs / 21.0, os / 16.0, ns / 7.0, ps / 3.0, ss / 1.0])

        # Don't need any other elements
        coa_completes = self._IsBalanced(atom_diff)
        # Stoichiometry right for CoA
        coa_completes &= numpy.all(n_coas == n_coas.astype(numpy.int))
        coa_completes &= numpy.all(n_coas ==
                                numpy.ones(n_coas.size) * n_coas[0])

        extra_coas = int(n_coas[0])
        if not coa_completes or not extra_coas:
            return None

        # Number of missing CoAs
        return extra_coas

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
        if from_id == to_id:
            return

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
            self.reactants[j].coeff += how_many
            self._Dedup()

        # clear the cache since the reaction has changed
        self._catalyzing_enzymes = None

    def _Dedup(self):
        """
            Collapses duplicate compounds and removes ones with
            coefficients of zero.
        """
        kegg_id_to_index = {}
        for i, c in enumerate(self.reactants):
            first_i = kegg_id_to_index.setdefault(c.compound.kegg_id, i)
            if i != first_i:
                self.reactants[first_i].coeff += c.coeff
                c.coeff = 0

        self.reactants = list(filter(lambda x: x.coeff != 0, self.reactants))

        # always make sure that H2O is the last reactant (so that it will
        # appear last in the chemical formula)
        i_h2o = self._FindCompoundIndex('C00001')
        if i_h2o is not None:
            self.reactants = self.reactants[:i_h2o] + \
                             self.reactants[(i_h2o + 1):] + \
                             [self.reactants[i_h2o]]

    def _FilterProtonsAndElectrons(self):
        """
            Removes Protons from the list of compounds.
            Since we use Bob Alberty's framework for biochemical reactions,
            there is no meaning for having H+ in a reaction.
        """
        self.reactants = list(filter(lambda c: c.compound.kegg_id not in
                                     ['C00080', 'C05359'], self.reactants))

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

    def TryBalanceWithCoA(self):
        """Try to balance the reaction with water.

        Returns:
            True if the reaction is balanced already or with
            additional waters on either side.
        """
        extra_coas = self._ExtraCoAs()
        if extra_coas is None:
            # cannot balance the reaction with CoA only
            return False
        if extra_coas != 0:
            self._AddCompound('C00010', extra_coas)
            self._Dedup()
        return True

    def TryBalanceWithPi(self):
        """Try to balance the reaction with water.

        TODO write a generic version of this?

        Returns:
            True if the reaction is balanced already or with
            additional waters on either side.
        """
        extra_pis = self._ExtraPis()
        if extra_pis is None:
            # cannot balance the reaction with Pi only
            return False
        if extra_pis != 0:
            self._AddCompound('C00009', extra_pis)
            self._Dedup()
        return True

    def CanBalanceWithWater(self):
        """Returns True if balanced with or without water."""
        try:
            extra_waters = self._ExtraWaters()
        except ReactantFormulaMissingError:
            return True

        return extra_waters is not None

    def CanBalanceWithCoA(self):
        """Returns True if balanced with or without CoA."""
        try:
            extra_coas = self._ExtraCoAs()
        except ReactantFormulaMissingError:
            return True

        return extra_coas is not None

    def CanBalanceWithPi(self):
        """Returns True if balanced with or without Pi."""
        try:
            extra_pis = self._ExtraPis()
        except ReactantFormulaMissingError:
            return True

        return extra_pis is not None

    def BalanceElectrons(self,
                         acceptor_id='C00003',           # NAD+
                         reduced_acceptor_id='C00004',
                         n_e=2):  # NADH
        """Try to balance the reaction electons.

        By default acceptor and reduced accepter differ by 2e-.
        """
        net_electrons = self._GetElectronDiff()
        if net_electrons != 0:
            self._AddCompound(reduced_acceptor_id, -net_electrons/float(n_e))
            self._AddCompound(acceptor_id, net_electrons/float(n_e))
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
        sum_logs = sum([c.coeff*numpy.log(c.phase.PhysiologicalValue())
                        for c in self.reactants])

        _r = constants.R
        _t = constants.DEFAULT_TEMP
        return _r * _t * sum_logs

    def DeltaG0Prime(self, aq_params=None):
        """Compute the DeltaG0' for a reaction.

        Args:
            aq_params: override whatever params are currently defined.

        Returns:
            The DeltaG0' for this reaction, or None if data was missing.
        """
        if self._dg0_prime is not None:
            logging.debug("Using cached dG0'")
            return self._dg0_prime

        logging.debug('Aqueous Params = ' + str(self.aq_params))
        aqp = aq_params or self.aq_params
        c_dg0_prime_list = [c.DeltaG0Prime(aqp) for c in self.reactants]

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
        if dg0_prime is None:
            return None
        return dg0_prime + correction

    def DeltaGPrime(self):
        """Compute the DeltaG' for a reaction.

        Returns:
            The DeltaG' for this reaction, or None if data was missing.
        """
        dg0_prime = self.DeltaG0Prime()
        correction = self._GetConcentrationCorrection()
        if dg0_prime is None:
            return None
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

        expo = numpy.floor(x)
        prefactor = 10**(x - expo)
        if abs(expo) <= 2:
            return '%.3g' % (10**x)
        else:
            return '%.1f &times; 10<sup>%d</sup>' % (prefactor, expo)

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
        if self._uncertainty is None:
            x, g = Preprocessing.GetReactionVectors(self.reactants)
            s_cc = Preprocessing.DeltaGUncertainty(x, g)
            logging.debug('s_cc = %g' % s_cc)
            self._uncertainty = 1.96*s_cc

            if Preprocessing.IsUsingGroupContributions(x, g):
                logging.debug('reaction is using GC')
                self.is_using_gc = True
            else:
                logging.debug('reaction is not using GC')
                self.is_using_gc = False

        return self._uncertainty

    def ExtraAtoms(self):
        try:
            diff = self._GetAtomDiff()
        except ReactantFormulaMissingError:
            return None
        diff.pop('H', 0)
        extras = list(filter(lambda t: t[1] > 0, diff.items()))
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
        short = list(filter(lambda t: t[1] < 0, diff.items()))
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

    def IsPhysiologicalConcentration(self):
        for c in self.reactants:
            if c.phase and not c.phase.IsPhysiological():
                return False
        return True

    def GetSourceReferenceLink(self, source_name):
        try:
            source = models.ValueSource.objects.get(name=source_name)
            url = source.url
        except Exception:
            url = "/data_refs"
        url = url or "/data_refs"
        return url

    def GetSourceReferences(self):
        """Returns a list of two-tuples (name, url).

        Assuming that all reactants in the chosen priority group have
        the same source reference, we provide the ref from the first
        reactant to represent all of them.
        """
        source_names = set(map(lambda x: str(x.compound._GetDGSource()),
                               self.reactants))
        urls = [self.GetSourceReferenceLink(n) for n in source_names]
        return zip(source_names, urls)

    def GetComponentContributionAnalysis(self):
        x, g = Preprocessing.GetReactionVectors(self.reactants)
        return Preprocessing.Analyze(x, g)

    substrates = property(GetSubstrates)
    products = property(GetProducts)
    is_reactant_formula_missing = property(IsReactantFormulaMissing)
    reactants_with_missing_formula = property(GetReactantFormulaMissing)
    is_empty = property(IsEmpty)
    is_balanced = property(IsBalanced)
    is_electron_balanced = property(IsElectronBalanced)
    balanced_with_water = property(CanBalanceWithWater)
    balanced_with_coa = property(CanBalanceWithCoA)
    balanced_with_pi = property(CanBalanceWithPi)
    extra_atoms = property(ExtraAtoms)
    missing_atoms = property(MissingAtoms)
    extra_electrons = property(ExtraElectrons)
    missing_electrons = property(MissingElectrons)
    special_reaction_warning = property(GetSpecialReactionWarning)
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
    link_url = property(GetHyperlink)
    ph_graph_link = property(GetPhGraphLink)
    pmg_graph_link = property(GetPMgGraphLink)
    is_graph_link = property(GetIonicStrengthGraphLink)
    catalyzing_enzymes = property(_GetCatalyzingEnzymes)
    is_phys_conc = property(IsPhysiologicalConcentration)
    source_references = property(GetSourceReferences)
    analyze_cc = property(GetComponentContributionAnalysis)


class StoredReaction(models.Model):
    """A reaction stored in the database."""
    # The ID of this reaction in KEGG.
    kegg_id = models.CharField(max_length=10, null=True)

    # a JSON representation of reaction (i.e. the reactants and coefficients)
    reactants = models.TextField(null=True)

    # a hash string for fast lookup of enzyme names by reaction
    reaction_hash = models.CharField(max_length=128, db_index=True)

    # a cache for the ToString() function which takes too long due to DB comm
    reaction_string = models.TextField(null=True)

    # a cache for the Link() function which takes too long due to DB comm
    link = models.TextField(null=True)

    @staticmethod
    def FromJson(rd):
        return StoredReaction(kegg_id=rd['RID'],
                              reactants=json.dumps(rd['reaction']))

    def GetSparseRepresentation(self):
        sparse = {}
        for coeff, kegg_id in json.loads(self.reactants):
            sparse[kegg_id] = coeff
        if 'C00080' in sparse:  # ignore H+ in stored reactions
            sparse.pop('C00080')
        return sparse

    @staticmethod
    def _CompoundToString(kegg_id, coeff):
        try:
            compound = apps.get_model('gibbs.Compound').objects.get(
                kegg_id=kegg_id)
            name = compound.FirstName()
        except Exception as e:
            logging.warning('Cannot find the name for %s' % kegg_id)
            logging.warning(str(e))
            name = kegg_id

        if coeff == 1:
            return name
        else:
            return "%g %s" % (coeff, name)

    def ToString(self):
        """
            String representation.
        """
        # TODO: need to replace the KEGG IDs with the common names of the
        #       compounds
        left = []
        right = []
        for coeff, kegg_id in json.loads(self.reactants):
            if coeff < 0:
                left.append(StoredReaction._CompoundToString(kegg_id, -coeff))
            elif coeff > 0:
                right.append(StoredReaction._CompoundToString(kegg_id, coeff))
        return "%s = %s" % (' + '.join(left), ' + '.join(right))

    @staticmethod
    def HashableReactionString(sparse):
        """Return a hashable string for a biochemical reaction.

        The string fully identifies the biochemical reaction up to
        directionality. If it is equal to another reaction's string,
        then they have identical stoichiometry up to their directionality.

        Args:
            sparse: a dictionary whose keys are kegg_id and values are
                    stoichiometric coefficients
        """
        if len(sparse) == 0:
            return ''

        # sort according to KEGG ID and normalize the stoichiometric
        # coefficients such that the coeff of the reactant with the lowest
        # ID will be 1
        kegg_id_list = sorted(sparse.keys())
        if sparse[kegg_id_list[0]] == 0:
            raise Exception('One of the stoichiometric coefficients is 0')
        norm_factor = 1.0 / sparse[kegg_id_list[0]]
        return ' + '.join(['%g %s' % (norm_factor*sparse[kegg_id], kegg_id)
                           for kegg_id in kegg_id_list])

    @staticmethod
    def HashReaction(sparse):
        md5 = hashlib.md5()
        md5.update(StoredReaction.HashableReactionString(sparse))
        return md5.hexdigest()

    @staticmethod
    def GetAtpHydrolysisHash():
        atp_sparse = {'C00002': -1, 'C00001': -1, 'C00008': 1, 'C00009': 1}
        return StoredReaction.HashableReactionString(atp_sparse)

    @staticmethod
    def GetCO2HydrationHash():
        co2_sparse = {'C00011': -1, 'C00001': -1, 'C00288': 1}
        return StoredReaction.HashableReactionString(co2_sparse)

    def GetHashableReactionString(self):
        return StoredReaction.HashableReactionString(
            self.GetSparseRepresentation())

    def GetHash(self):
        return StoredReaction.HashReaction(self.GetSparseRepresentation())

    def GenerateHash(self):
        self.reaction_hash = self.GetHash()
        self.reaction_string = self.ToString()
        self.link = self.Link()

    def __str__(self):
        """String representation."""
        return self.ToString()

    def Link(self):
        """
            Returns a link to this reaction's page.
        """
        try:
            rxn = self.ToReaction()
            return rxn.GetHyperlink(self.ToString())
        except AttributeError:
            raise Exception('Cannot find one of the compounds in the database')

    def ToReaction(self, priority=1, aq_params=None):
        """
            returns a list of CSV rows with the following columns:
            kegg ID, dG0_prime, pH, ionic_strength, T, Note
        """
        reactants = [CompoundWithCoeff.FromId(coeff, kegg_id)
                     for coeff, kegg_id in json.loads(self.reactants)]
        rxn = Reaction(reactants)
        return rxn


class Enzyme(models.Model):
    """A single enzyme."""
    # EC class enzyme.
    ec = models.CharField(max_length=10)

    # A list of common names of the enzyme, used for searching.
    common_names = models.ManyToManyField(CommonName)

    # List of reactions this enzyme catalyzes.
    reactions = models.ManyToManyField(StoredReaction)

    def HasData(self):
        """Checks if it has enough data to display."""
        return self.ec and self.reactions.all()

    def ToJson(self):
        """Returns as a JSON-friendly object."""
        return {'EC': self.ec,
                'name': str(self.first_name)}

    def Link(self):
        """Returns a link to this reactions page."""
        return '/enzyme?ec=%s' % self.ec

    def KeggLink(self):
        """Returns a link to the KEGG page for this enzyme."""
        return 'http://kegg.jp/dbget-bin/www_bget?ec:%s' % self.ec

    def BrendaLink(self):
        """Returns a link to the BRENDA page for this enzyme."""
        return 'http://www.brenda-enzymes.org/php/result_flat.php4?ecno=%s' \
            % self.ec

    def AllReactions(self):
        """Returns all the reactions."""
        return self.reactions.all()

    def FirstName(self):
        """The first name in the list of names."""
        return self.all_common_names[0]

    def __unicode__(self):
        """Return a single string identifier of this enzyme."""
        return str(self.FirstName())

    all_common_names = property(lambda self: self.common_names.all())
    all_cofactors = property(lambda self: self.cofactors.all())
    first_name = property(FirstName)
    all_reactions = property(AllReactions)
    kegg_link = property(KeggLink)
    brenda_link = property(BrendaLink)
    link = property(Link)
