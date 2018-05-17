from util.SBtab.SBtabDict import SBtabDict
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from io import StringIO
import logging
from util import constants
from gibbs import service_config
from django.apps import apps
from pathway.bounds import Bounds
from os import path
from equilibrator.settings import BASE_DIR
from gibbs.conditions import AqueousParams

class PathwayParseError(Exception):
    
    def __init__(self, message):
        super(PathwayParseError, self).__init__(message)
        self.message = message

class InvalidReactionFormula(PathwayParseError):
    pass


class UnbalancedReaction(PathwayParseError):
    pass


class ViolatesFirstLaw(PathwayParseError):
    pass

def HtmlConcentration(conc):
    if conc <= 9.999e-4:
        return '%.1f &mu;M' % (1e6*conc)
    if conc <= 0.9999:
        return '%.1f mM' % (1e3*conc)
    return '%.1f M' % conc


# TODO - implement all the functions common to MDF and ECM here, such as
# reading the reaction list, storing the dG'0 values, the pH and ionic strength,
# concentration bounds, etc.
class ParsedPathway(object):
    
    EXPECTED_TNAMES = []
    COFACTORS_FNAME = path.join(BASE_DIR, 'pathway/data/cofactors.csv')
    DEFAULT_BOUNDS = Bounds.from_csv_filename(
        COFACTORS_FNAME, default_lb=1e-6, default_ub=0.1)
    
    def __init__(self, reactions, fluxes, bounds=None, aq_params=None):
        """Initialize.

        Args:
            reactions: a list of gibbs.reaction.Reaction objects.
            fluxes: np.array of relative fluxes in same order as reactions.
            bounds: bounds on metabolite concentrations.
                Uses default bounds if None provided.
            aq_params: specify the pH, ionic strength, etc. at which the
                dG values are calculated. May be omitted.
        """
        assert len(reactions) == len(fluxes)

        self._reactions = reactions
        self.reaction_kegg_ids = [r.stored_reaction_id for r in reactions]
        self.aq_params = AqueousParams()  # Default values

        self._fluxes = np.array(fluxes)

        self.bounds = bounds or ParsedPathway.DEFAULT_BOUNDS

        self.S, self.compound_kegg_ids = self._build_stoichiometric_matrix()
        self.compounds_by_kegg_id = self._get_compounds()
        self._compounds = [self.compounds_by_kegg_id[cid]
                          for cid in self.compound_kegg_ids]

        nr, nc = self.S.shape

        net_rxn_stoich = (self._fluxes.reshape((nr, 1)) * self.S).sum(axis=0)
        net_rxn_data = []
        for coeff, kid in zip(net_rxn_stoich, self.compound_kegg_ids):
            if coeff != 0:
                net_rxn_data.append(self._reactant_dict(coeff, kid))
        self.net_reaction = apps.get_model('gibbs.reaction').FromIds(net_rxn_data, fetch_db_names=True)
    
    @classmethod
    def from_sbtab_file(cls, f):
        sbtabs = SBtabDict.FromSBtabFile(f)
        if not set(cls.EXPECTED_TNAMES).issubset(sbtabs.keys()):
            raise PathwayParseError('Make sure the pathway model SBtab file '
                                    'contains all necessary tables: ' + 
                                    ', '.join(cls.EXPECTED_TNAMES))
        return cls.from_sbtabs(sbtabs)

    def is_empty(self):
        return len(self.reactions) == 0
        
    def analyze(self):
        raise Exception('ParsedPathway is a virtual class, one should not use its '
                        'static functions.')

    @property
    def reactions(self):
        return self._reactions

    @property
    def compounds(self):
        return self._compounds

    @property
    def fluxes(self):
        return self._fluxes

    @staticmethod
    def _reactant_dict(coeff, kid, negate=False):
        """Returns dictionary format expected by Reaction.FromIds."""
        if negate:
            coeff = -1*coeff
        d = {'kegg_id': kid, 'coeff': coeff, 'name': kid,
             'phase': constants.AQUEOUS_PHASE_NAME}
        if kid == 'C00001':
            # Water is not aqueous. Hate that this is hardcoded.
            d['phase'] = constants.LIQUID_PHASE_NAME
        return d

    @classmethod
    def get_data_from_sbtabs(cls, sbtabs):
        query_parser = service_config.Get().query_parser
        reaction_df = sbtabs['Reaction'].toDataFrame()
        bounds_df = sbtabs['ConcentrationConstraint'].toDataFrame()
        
        name_to_cid = dict(
            zip(bounds_df['Compound'],
                bounds_df['Compound:Identifiers:kegg.compound']))

        reactions = []
        for idx in reaction_df.index:
            row = reaction_df.loc[idx]
            rxn_formula = row['ReactionFormula']
            logging.debug(rxn_formula)
            parsed_rxn = query_parser.ParseReactionQuery(rxn_formula)

            rxn_ds = []
            for coeff, name in parsed_rxn.substrates:
                cid = name_to_cid[name]
                rxn_ds.append(cls._reactant_dict(coeff, cid, negate=True))
            for coeff, name in parsed_rxn.products:
                cid = name_to_cid[name]
                rxn_ds.append(cls._reactant_dict(coeff, cid, negate=False))
            rxn = apps.get_model('gibbs.reaction').FromIds(rxn_ds, fetch_db_names=True)

            if not rxn.IsBalanced():
                raise UnbalancedReaction(
                    "ReactionFormula '%s' is not balanced" % rxn_formula)
            if not rxn.IsElectronBalanced():
                raise UnbalancedReaction(
                    "ReactionFormula '%s' is not redox balanced" % rxn_formula)

            reactions.append(rxn)

        reaction_ids = reaction_df['ID']
        
        # read flux table
        if 'RelativeFlux' in sbtabs:
            flux_df = sbtabs['RelativeFlux'].toDataFrame()
        elif 'Flux' in sbtabs:
            flux_df = sbtabs['Flux'].toDataFrame()
        fluxes = flux_df[flux_df['QuantityType'] == 'flux']
        reaction_fluxes = dict(zip(fluxes['Reaction'], fluxes['Value']))
        fluxes_ordered = [float(reaction_fluxes[rid]) for rid in reaction_ids]

        bounds = Bounds.from_sbtab(sbtabs['ConcentrationConstraint'])
        
        return reactions, fluxes_ordered, bounds

    @classmethod
    def from_sbtabs(cls, sbtabs):
        raise Exception('ParsedPathway is a virtual class, one should not use its '
                        'static functions.')

    def _get_compounds(self):
        """Returns a dictionary of compounds by KEGG ID."""
        compounds = {}
        for r in self.reactions:
            for cw_coeff in r.reactants:
                c = cw_coeff.compound
                compounds[c.kegg_id] = c
        return compounds

    def _build_stoichiometric_matrix(self):
        """Builds a stoichiometric matrix.

        Returns:
            Two tuple (S, compounds) where compounds is the KEGG IDs of the compounds
            in the order defining the column order of the stoichiometric matrix S.
        """
        compounds = []
        sparses = []
        for r in self.reactions:
            s = r.GetSparseRepresentation()
            sparses.append(s)
            for kegg_id in s:
                compounds.append(kegg_id)
        compounds = sorted(set(compounds))

        # reactions on the rows, compounds on the columns
        n_reactions = len(self.reactions)
        n_compounds = len(compounds)
        smat = np.zeros((n_reactions, n_compounds))
        for i, s in enumerate(sparses):
            for j, c in enumerate(compounds):
                smat[i, j] = s.get(c, 0)
        return smat, compounds

    
class ReactionData(object):

    def __init__(self, reaction, flux, dGr=0, shadow_price=0, enz_conc=0):
        """
        Args:
            reaction: kegg reaction object.
                should be set to contain user-defined dG0
            flux: amount of relative flux in pathway.
            dGr: dG in MDF conditions.
            shadow_price: shadow price associated with this rxn.
        """
        self.reaction = reaction
        self.flux = flux
        self.dGr = dGr
        self.shadow_price = shadow_price
        self.enz_conc = enz_conc

    @property
    def dG0_prime(self):
        return self.reaction.dg0_prime

    @property
    def dGm_prime(self):
        return self.reaction.dgm_prime
    
    @property
    def html_enzyme_concentration(self):
        return HtmlConcentration(self.enz_conc)

class CompoundData(object):
    def __init__(self, compound, concentration_bounds,
                 concentration=0, shadow_price=0, cost=0):
        self.compound = compound
        self.concentration = concentration
        self.shadow_price = shadow_price
        self.lb, self.ub = concentration_bounds
        self.cost = cost

    @property
    def compound_name(self):
        return self.compound.name.name

    @property
    def is_water(self):
        return self.compound.kegg_id == 'C00001'

    @property
    def link_url(self):
        return '/compound?compoundId=%s' % self.compound.kegg_id

    @property
    def bounds_equal(self):
        return self.lb == self.ub

    @property
    def html_concentration(self):
        return HtmlConcentration(self.concentration)

    @property
    def html_lb(self):
        return HtmlConcentration(self.lb)

    @property
    def html_ub(self):
        return HtmlConcentration(self.ub)

class PathwayAnalysisData(object):
    
    def __init__(self, parsed_pathway):
        self.parsed_pathway = parsed_pathway
        
        rxns = parsed_pathway.reactions
        fluxes = parsed_pathway.fluxes
        self.reaction_data = [
            ReactionData(*t) for t in zip(rxns, fluxes)]

        compounds = parsed_pathway.compounds
        cbounds = map(parsed_pathway.bounds.GetBoundTuple,
                      parsed_pathway.compound_kegg_ids)
        
        self.compound_data = [
            CompoundData(*t) for t in zip(compounds, cbounds)]

    @property
    def score(self):
        return 0
    
    @property
    def is_mdf(self):
        return False
    
    @property
    def is_ecm(self):
        return False

    @property
    def reaction_plot_svg(self):
        pass

    @property
    def metabolite_plot_svg(self):
        default_lb = self.parsed_pathway.bounds.default_lb
        default_ub = self.parsed_pathway.bounds.default_ub

        concs = []
        lbs = []
        ubs = []
        cnames = []
        shadow_prices = []
        
        for c in reversed(self.compound_data):
            if c.is_water:
                continue
            concs.append(c.concentration)
            lbs.append(c.lb)
            ubs.append(c.ub)
            cnames.append(c.compound_name)
            shadow_prices.append(c.shadow_price)
        
        Nc = len(cnames)
        
        concs = np.array(concs)
        lbs = np.array(lbs)
        ubs = np.array(ubs)
        shadow_prices = np.array(shadow_prices)
        ys = np.arange(0, Nc)

        bounds_equal = np.where(lbs == ubs)
        ys_equal = ys[bounds_equal]
        concs_equal = concs[bounds_equal]

        conc_figure = plt.figure(figsize=(8, 6))
        sns.set_style('darkgrid')
        plt.axes([0.2, 0.1, 0.9, 0.9])
        plt.axvspan(1e-8, default_lb, color='y', alpha=0.5)
        plt.axvspan(default_ub, 1e3, color='y', alpha=0.5)
        plt.scatter(concs, ys, figure=conc_figure,
                    label='Variable Concentrations')
        plt.scatter(concs_equal, ys_equal, figure=conc_figure, color='y',
                    label='Fixed Concentrations')
        
        # Special color for metabolites with nonzero shadow prices.
        nz_shadow = np.where(shadow_prices != 0)
        ys_nz_shadow = ys[nz_shadow]
        concs_nz_shadow = concs[nz_shadow]
        plt.scatter(concs_nz_shadow, ys_nz_shadow, figure=conc_figure,
                    color='r', label='Variable Concentrations')

        plt.xticks(family='sans-serif', figure=conc_figure)
        plt.yticks(ys, cnames, family='sans-serif',
            fontsize=8, figure=conc_figure)
        plt.xlabel('Concentration (M)', family='sans-serif',
            figure=conc_figure)
        plt.xscale('log')

        plt.xlim(1e-7, 1.0)
        plt.ylim(-0.5, Nc - 0.5)

        svg_data = StringIO()
        conc_figure.savefig(svg_data, format='svg')
        return svg_data.getvalue()
