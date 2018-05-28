from util.SBtab.SBtabDict import SBtabDict, SBtabError
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
import pandas as pd
import csv
from pathway.concs import ConcentrationConverter
from scipy import linalg

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
    if conc <= 9.999e-7:
        return '%.1f nM' % (1e9*conc)
    if conc <= 9.999e-4:
        return '%.1f &mu;M' % (1e6*conc)
    if conc <= 0.9999:
        return '%.1f mM' % (1e3*conc)
    return '%.1f M' % conc


class ParsedPathway(object):

    EXPECTED_TNAMES = set(['Reaction', 'Compound', 'Parameter', 'Flux',
                           'ConcentrationConstraint'])
    COFACTORS_FNAME = path.join(BASE_DIR, 'pathway/data/cofactors.csv')
    DEFAULT_BOUNDS = Bounds.from_csv_filename(
        COFACTORS_FNAME, default_lb=1e-6, default_ub=0.1)
    
    def __init__(self, reactions, fluxes, bounds=None, aq_params=None,
                 reaction_ids=None):
        """Initialize.

        Args:
            reactions: a list of gibbs.reaction.Reaction objects.
            fluxes: np.array of relative fluxes in same order as reactions.
            bounds: bounds on metabolite concentrations.
                Uses default bounds if None provided.
            aq_params: specify the pH, ionic strength, etc. at which the
                dG values are calculated. May be omitted.
            reaction_ids: optional way to provide names for reactions that 
                will be used in the report later (otherwise, short names will
                be generated automatically)
        """
        assert len(reactions) == len(fluxes)

        self.reactions = reactions
        if reaction_ids is None:
            self.reaction_ids = []
            for i, rxn in enumerate(self.reactions):
                self.reaction_ids.append(rxn.GenerateUniqueName(i))
        else:
            self.reaction_ids = reaction_ids

        self.aq_params = aq_params or AqueousParams()  # Default values

        self.fluxes = np.array(fluxes)

        self.bounds = bounds or ParsedPathway.DEFAULT_BOUNDS

        self.S, self.compound_kegg_ids = self._build_stoichiometric_matrix()
        self.compounds_by_kegg_id = self.get_compounds()
        self.compounds = [self.compounds_by_kegg_id[cid]
                          for cid in self.compound_kegg_ids]

        nr, nc = self.S.shape

        net_rxn_stoich = (self.fluxes.reshape((nr, 1)) * self.S).sum(axis=0)
        net_rxn_data = []
        for coeff, kid in zip(net_rxn_stoich, self.compound_kegg_ids):
            if coeff != 0:
                net_rxn_data.append(ParsedPathway.reactant_dict(coeff, kid))
        self.net_reaction = apps.get_model('gibbs.reaction').FromIds(
                net_rxn_data, fetch_db_names=True)

    @staticmethod
    def reactant_dict(coeff, kid, negate=False):
        """Returns dictionary format expected by Reaction.FromIds."""
        if negate:
            coeff = -1*coeff
        d = {'kegg_id': kid, 'coeff': coeff, 'name': kid,
             'phase': constants.AQUEOUS_PHASE_NAME}
        if kid == 'C00001':
            # Water is not aqueous. Hate that this is hardcoded.
            d['phase'] = constants.LIQUID_PHASE_NAME
        return d

    def validate_dGs(self):
        """
            dGr should be orthogonal to nullspace of S
            If not, dGr is not contained in image(S) and then there
            is no consistent set of dGfs that generates dGr and the
            first law of thermo is violated by the model.
        """
        Smat = np.matrix(self.S)
        Spinv = linalg.pinv(Smat)
        null_proj = np.matrix(np.eye(Smat.shape[0])) - Smat*Spinv
        projected = null_proj * self.dG0_r_primes
        return np.all(projected < 1e-8)

    def is_empty(self):
        return len(self.reactions) == 0
        
    def analyze(self):
        raise Exception('ParsedPathway is a virtual class, one should not use its '
                        'static functions.')

    @property
    def dG0_r_primes(self):
        return np.matrix([rxn.dg0_r_prime for rxn in self.reactions]).T

    @property
    def reactions_balanced(self):
        """Returns true if all pathway reactions are electron and atom-wise balanced."""
        atom_balaned = [r.IsBalanced() for r in self.reactions]
        electron_balaned = [r.IsElectronBalanced() for r in self.reactions]

        balanced = np.logical_and(atom_balaned, electron_balaned)
        return np.all(balanced)

    def print_reactions(self):
        for f, r in zip(self.fluxes, self.reactions):
            print('%sx %s' % (f, r))

    @staticmethod
    def from_sbtab(sbtabs):
        """
            Initializes a ParsedPathway object with data a SBtabDict
        """
        name_to_cid = sbtabs.GetDictFromTable('Compound', 'ID',
                                              'Identifiers:kegg.compound')
        query_parser = service_config.Get().query_parser
            
        def formula2rxn(formula):
            parsed_rxn = query_parser.ParseReactionQuery(formula)

            rxn_ds = []
            for coeff, name in parsed_rxn.substrates:
                cid = name_to_cid[name]
                rxn_ds.append(ParsedPathway.reactant_dict(coeff, cid, negate=True))
            for coeff, name in parsed_rxn.products:
                cid = name_to_cid[name]
                rxn_ds.append(ParsedPathway.reactant_dict(coeff, cid, negate=False))
            rxn = apps.get_model('gibbs.reaction').FromIds(rxn_ds, fetch_db_names=True)

            if not rxn.IsBalanced():
                raise UnbalancedReaction(
                    "ReactionFormula '%s' is not balanced" % formula)
            if not rxn.IsElectronBalanced():
                raise UnbalancedReaction(
                    "ReactionFormula '%s' is not redox balanced" % formula)
            return rxn
        
        # read an parse reactions, creating equilibrator Reaction objects
        reactions = list(map(formula2rxn,
                             sbtabs.GetColumnFromTable('Reaction',
                                                       'ReactionFormula')))
        reaction_ids = sbtabs.GetColumnFromTable('Reaction', 'ID')
        
        # read flux table and sort it according to the order of reactions
        # in the reaction table
        reaction_fluxes = sbtabs.GetDictFromTable('Flux', 'Reaction', 'Value')
        fluxes = [float(x) for x in map(reaction_fluxes.get, reaction_ids)]

        bounds = Bounds.from_sbtab(sbtabs['ConcentrationConstraint'])
        
        # read the general aqueous parameters from the ReactionConstant
        # table header
        parameter_sbtab = sbtabs['Parameter']
        aq_params = AqueousParams()  # Default values
        try:
            pH, ionic_strength, ionic_strength_units = \
                map(parameter_sbtab.getCustomTableInformation,
                    ['pH', 'IonicStrength', 'IonicStrengthUnit'])
            aq_params.pH = float(pH)
            c = float(ionic_strength)
            c = ConcentrationConverter.to_molar_string(c, ionic_strength_units)
            aq_params.ionic_strength = c
        except SBtabError:
            logging.debug('pH or I unspecified in SBtab, using default values')
        
        # get the equilibrium constants from the Parameters table
        keqs_df = parameter_sbtab.toDataFrame()
        keqs = keqs_df[keqs_df.QuantityType == 'equilibrium constant'].set_index('Reaction')
        reaction_keqs = np.array(list(map(float, keqs.loc[reaction_ids, 'Value'])))
        dG0_r_primes = -constants.RT * np.log(reaction_keqs)

        # Manually set the delta G values on the reaction objects
        for rxn, dg in zip(reactions, dG0_r_primes):
            rxn.dg0_r_prime = dg

        return ParsedPathway(reactions, fluxes, bounds, aq_params, reaction_ids)

    def get_compounds(self):
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
    
    @classmethod
    def from_csv(cls, f, bounds=None, aq_params=None):
        """
            This is the superclass implementation, that should be used
            by the inherited classes at the beginning of from_csv()
        """
        rxn_matcher = service_config.Get().reaction_matcher
        query_parser = service_config.Get().query_parser
        aq_params = aq_params or AqueousParams()

        reaction_df = pd.read_csv(f,
                                  dtype={'ReactionFormula':str, 'Flux':float})
        
        if len(reaction_df.columns) != 2:
            raise InvalidReactionFormula(
                "Input CSV file must have exactly 2 columns")
        if reaction_df.columns[0] != 'ReactionFormula':
            raise InvalidReactionFormula(
                "First column in CSV file must be 'ReactionFormula'")
        if reaction_df.columns[1] != 'Flux':
            raise InvalidReactionFormula(
                "Second column in CSV file must be 'Flux'")
        
        fluxes = reaction_df.Flux.fillna(0.0).tolist()
        
        reactions = []
        for formula in reaction_df.ReactionFormula:
            if not formula:
                raise InvalidReactionFormula('Found empty ReactionFormula')

            logging.debug('formula = %f x (%s)', formula)

            if not query_parser.IsReactionQuery(formula):
                raise InvalidReactionFormula("Failed to parse '%s'", formula)

            parsed = query_parser.ParseReactionQuery(formula)

            matches = rxn_matcher.MatchReaction(parsed)
            best_match = matches.GetBestMatch()
            rxn = apps.get_model('gibbs.reaction').FromIds(
                best_match, fetch_db_names=True)

            if not rxn.IsBalanced():
                raise UnbalancedReaction(
                    "ReactionFormula '%s' is not balanced" % formula)
            if not rxn.IsElectronBalanced():
                raise UnbalancedReaction(
                    "ReactionFormula '%s' is not redox balanced" % formula)

            # calculate the dG'0 value for this reaction and store it in 
            # a float variable (so it will be fast to access it later)
            rxn.dg0_r_prime = rxn.DeltaG0Prime(aq_params)
            reactions.append(rxn)

        return ParsedPathway(reactions, fluxes, bounds, aq_params)
    
    @staticmethod
    def create_sbtab_header(tablename, tabletype, **kwargs):
        s = "!!SBtab TableName='%s' TableType='%s'" % (tablename, tabletype)
        s += " Document='Pathway Model' SBtabVersion='1.0'"
        for name, value in kwargs.items():
            s += " %s='%s'" % (name, value)
        s += "\n"
        return s
            
    def to_sbtab(self):
        """
            Returns a full SBtab description of the model.

            Description includes reaction fluxes and per-compound bounds.
            
            This is the superclass implementation, that should be used
            by the inherited classes at the beginning of to_sbtab()
        """
        # Reaction table
        reaction_header = ParsedPathway.create_sbtab_header('Reaction', 'Reaction')
        reaction_cols = ['!ID', '!Name', '!ReactionFormula', '!Identifiers:kegg.reaction']

        sio = StringIO()
        sio.write(reaction_header)
        writer = csv.DictWriter(sio, reaction_cols, dialect='excel-tab')
        writer.writeheader()

        for rxn_id, rxn in zip(self.reaction_ids, self.reactions):
            d = {'!ID': rxn_id,
                 '!Name': rxn_id,
                 '!ReactionFormula': rxn.GetSlugQueryString(),
                 '!Identifiers:kegg.reaction': rxn.stored_reaction_id}
            writer.writerow(d)
        sio.write('%\n')
        
        # Compound table
        compound_header = ParsedPathway.create_sbtab_header('Compound', 'Compound')
        compound_cols = ['!ID', '!Name', '!Identifiers:kegg.compound', '!IsConstant']
        sio.write(compound_header)
        writer = csv.DictWriter(sio, compound_cols, dialect='excel-tab')
        writer.writeheader()
        for cid, compound in self.compounds_by_kegg_id.items():
            d = {'!ID': compound.name_slug,
                 '!Name': compound.name.name,
                 '!Identifiers:kegg.compound': cid,
                 '!IsConstant': 'False'}
            writer.writerow(d)
        sio.write('%\n')

        # Flux table
        flux_header = ParsedPathway.create_sbtab_header('Flux', 'Quantity', Unit='mM/s')
        flux_cols = ['!QuantityType', '!Reaction', '!Reaction:Identifiers:kegg.reaction', '!Value']
        sio.write(flux_header)
        writer = csv.DictWriter(sio, flux_cols, dialect='excel-tab')
        writer.writeheader()

        for rxn_id, rxn, flux in zip(self.reaction_ids, self.reactions, self.fluxes):
            d = {'!QuantityType': 'flux',
                 '!Reaction': rxn_id,
                 '!Reaction:Identifiers:kegg.reaction': rxn.stored_reaction_id,
                 '!Value': flux}
            writer.writerow(d)
        sio.write('%\n')

        # ConcentrationConstraint table
        conc_header = ParsedPathway.create_sbtab_header('ConcentrationConstraint',
                                                        'Quantity', Unit='M')
        conc_cols = ['!QuantityType', '!Compound',
                     '!Compound:Identifiers:kegg.compound',
                     '!Concentration:Min', '!Concentration:Max']
        sio.write(conc_header)

        writer = csv.DictWriter(sio, conc_cols, dialect='excel-tab')
        writer.writeheader()
        for cid, compound in self.compounds_by_kegg_id.items():
            d = {'!QuantityType': 'concentration',
                 '!Compound': str(compound.name_slug),
                 '!Compound:Identifiers:kegg.compound': cid,
                 '!Concentration:Min': self.bounds.GetLowerBound(cid),
                 '!Concentration:Max': self.bounds.GetUpperBound(cid)}
            writer.writerow(d)
        sio.write('%\n')

        return sio.getvalue()
    
class PathwayAnalyzer(object):
    
    """
        A class for handling whole pathways, including the stoichiometry,
        naming conventions, and several standard parameters that are typically
        used in pathway analysis: concentration bounds, fluxes, aqueous 
        enviroment parameters, etc.
    """

    def __init__(self, parsed_pathway):
        self._parsed_pathway = parsed_pathway

    @property
    def reactions(self):
        return self._parsed_pathway.reactions

    @property
    def compounds(self):
        return self._parsed_pathway.compounds

    @property
    def fluxes(self):
        return self._parsed_pathway.fluxes

    @property
    def aq_params(self):
        return self._parsed_pathway.aq_params

    @property
    def reaction_ids(self):
        return self._parsed_pathway.reaction_ids

    @property
    def net_reaction(self):
        return self._parsed_pathway.net_reaction

    def to_sbtab(self):
        return self.parsed_pathway.to_sbtab()

    @classmethod
    def validate_sbtab(cls, sbtab):
        missing = ParsedPathway.EXPECTED_TNAMES.difference(sbtab.keys())
        if missing:
            raise PathwayParseError('Make sure the pathway model SBtab file '
                                    'contains these tables: ' + 
                                    ', '.join(missing))
        
    @classmethod
    def from_sbtab(cls, sbtabs):
        return PathwayAnalyzer(ParsedPathway.from_sbtab(sbtabs))

    @classmethod
    def from_sbtab_file(cls, fp):
        sbtab = SBtabDict.FromSBtabFile(fp)
        cls.validate_sbtab(sbtab)
        return cls.from_sbtab(sbtab)

    def validate_dGs(self):
        return self._parsed_pathway.validate_dGs()

    def is_empty(self):
        return self._parsed_pathway.is_empty()
    
    
class ReactionData(object):
    """
        A class for storing reaction-related results from pathway analysis,
        such as driving force, efficiency, etc.
    """

    def __init__(self, reaction, flux, name,
                 dGr=0, shadow_price=0, enz_conc=0,
                 min_enz_conc=0, eta_th=0, eta_sat=0):
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
        self.name = name
        self.dGr = dGr
        self.shadow_price = shadow_price
        self.enz_conc = enz_conc
        self.min_enz_conc = min_enz_conc
        self.eta_th = eta_th
        self.eta_sat = eta_sat

    @property
    def reaction_formula(self):
        return self.reaction.GetQueryString()

    @property
    def dG0_prime(self):
        return self.reaction.dg0_prime

    @property
    def dGm_prime(self):
        return self.reaction.dgm_prime
    
    @property
    def html_enzyme_concentration(self):
        return HtmlConcentration(self.enz_conc)
    
    @property
    def html_enzyme_minimum_conc(self):
        return HtmlConcentration(self.min_enz_conc)

class CompoundData(object):
    """
        A class for storing compound-related results from pathway analysis,
        such as concentrations, bounds, etc.
    """

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
    """
        A class for storing all the results of a pathway analysis.
        
    """
    
    COLOR_FIXED_CONCENTRATION = sns.color_palette('muted')[1]
    COLOR_VARIABLE_CONCENTRATION = sns.color_palette('muted')[0]
    COLOR_BOTTLENECK_CONCENTRATION = sns.color_palette('muted')[2]
    COLOR_GENERAL_CONCENTRATION = sns.color_palette('muted')[1]
    COLOR_GRID = '#999999'
    COLOR_DELTA_G_M = '#999999'
    COLOR_DELTA_G_MDF = sns.color_palette('muted')[0]
    COLOR_BOTTLENECK_REACTIONS = sns.color_palette('muted')[2]

    def __init__(self, parsed_pathway):
        self.parsed_pathway = parsed_pathway
        
        rxns = parsed_pathway.reactions
        fluxes = parsed_pathway.fluxes
        rids = parsed_pathway.reaction_ids
        self.reaction_data = [
            ReactionData(*t) for t in zip(rxns, fluxes, rids)]

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
    def reaction_names(self):
        return list(map(lambda rxn: rxn.name, self.reaction_data))

    @property
    def metabolite_plot_svg(self):
        """
            Generate the optimized concentrations graph, where compounds
            are ordered along the y-axis and their concentration is
            indicated by a point along the x-axis (in log-scale).
            The colors indicate whether this metabolite has a fixed
            concentration (i.e. lower bound = upper bound) and also
            if it has a non-zero shadow price (only in MDF).
        """
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

        with sns.axes_style('darkgrid'):
            conc_fig, ax = plt.subplots(1, 1, figsize=(8, 6))
            ax.grid(color=self.COLOR_GRID, linestyle='--', linewidth=1, alpha=0.2)

            ax.axvspan(1e-8, default_lb, color=self.COLOR_GENERAL_CONCENTRATION, alpha=0.5)
            ax.axvspan(default_ub, 1e3, color=self.COLOR_GENERAL_CONCENTRATION, alpha=0.5)
            ax.scatter(concs, ys, 
                       color=self.COLOR_VARIABLE_CONCENTRATION,
                       label='Variable concentrations')
            ax.scatter(concs_equal, ys_equal,
                       color=self.COLOR_FIXED_CONCENTRATION,
                       label='Fixed concentrations')
            
            # Special color for metabolites with nonzero shadow prices.
            nz_shadow = np.where(shadow_prices != 0)
            ys_nz_shadow = ys[nz_shadow]
            concs_nz_shadow = concs[nz_shadow]
            ax.scatter(concs_nz_shadow, ys_nz_shadow,
                       color=self.COLOR_BOTTLENECK_CONCENTRATION,
                       label='Bottleneck concentrations')
    
            ax.set_yticks(ys)
            ax.set_yticklabels(cnames)
            ax.set_xlabel('Concentration (M)')
            ax.set_xscale('log')
   
            ax.set_xlim(1e-7, 1.0)
            ax.set_ylim(-0.5, Nc - 0.5)
            ax.legend(loc='best', framealpha=0.5)
            conc_fig.tight_layout()
    
            svg_data = StringIO()
            conc_fig.savefig(svg_data, format='svg')
        return svg_data.getvalue()
