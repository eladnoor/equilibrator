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
import pandas as pd
import csv
from django.utils.text import slugify
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
    SBTAB_GENERIC_HEADER = \
        "!!SBtab TableName='%s' TableType='%s' Document='Pathway Model' SBtabVersion='1.0'"
    
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
        self.reaction_ids = []
        for i, rxn in enumerate(self.reactions):
            kegg_id = rxn.stored_reaction_id
            if rxn.catalyzing_enzymes:
                enz = str(rxn.catalyzing_enzymes[0].FirstName().name)
                enz_slug = slugify(enz)[:10]
                enz_slug = enz_slug.replace('-', '_')
                rxn_id = '%s_%s' % (enz_slug, kegg_id)
            elif not kegg_id:
                rxn_id = 'RXN%03d' % i
            self.reaction_ids.append(rxn_id)

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

    @classmethod
    def from_sbtab_file(cls, f):
        sbtabs = SBtabDict.FromSBtabFile(f)
        if not set(cls.EXPECTED_TNAMES).issubset(sbtabs.keys()):
            raise PathwayParseError('Make sure the pathway model SBtab file '
                                    'contains all necessary tables: ' + 
                                    ', '.join(cls.EXPECTED_TNAMES))
        return cls.from_sbtab(sbtabs)

    def is_empty(self):
        return len(self.reactions) == 0
        
    def analyze(self):
        raise Exception('ParsedPathway is a virtual class, one should not use its '
                        'static functions.')

    @property
    def reactions(self):
        return self._reactions

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
    def _from_sbtab(cls, sbtabs):
        """
            This is the superclass implementation, that should be used
            by the inherited classes at the beginning of from_sbtab()
        """
        name_to_cid = sbtabs.GetDictFromTable('Compound', 'ID',
                                              'Identifiers:kegg.compound')
        query_parser = service_config.Get().query_parser
            
        def formula2rxn(formula):
            parsed_rxn = query_parser.ParseReactionQuery(formula)

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
        
        reaction_ids = sbtabs['Reaction'].toDataFrame()['ID']

        # read the general aqueous parameters from the ReactionConstant
        # table header
        parameter_sbtab = sbtabs['Parameter']
        pH, ionic_strength, ionic_strength_units = \
            map(parameter_sbtab.getCustomTableInformation,
                ['pH', 'IonicStrength', 'IonicStrengthUnit'])
        aq_params = AqueousParams()  # Default values
        if pH:
            aq_params.pH = float(pH)
        if ionic_strength:
            c = float(ionic_strength)
            c = ConcentrationConverter.to_molar_string(c, ionic_strength_units)
            aq_params.ionic_strength = c

        # get the equilibrium constants from the Parameters table
        keqs_df = parameter_sbtab.toDataFrame()
        keqs = keqs_df[keqs_df.QuantityType == 'equilibrium constant'].set_index('Reaction')
        reaction_keqs = np.array(list(map(float, keqs.loc[reaction_ids, 'Value'])))
        dG0_r_primes = -constants.RT * np.log(reaction_keqs)

        # Manually set the delta G values on the reaction objects
        for rxn, dg in zip(reactions, dG0_r_primes):
            rxn.dg0_r_prime = dg

        return reactions, fluxes, bounds, aq_params

    @classmethod
    def from_sbtab(cls, sbtabs):
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
    
    @classmethod
    def _from_csv(cls, f, bounds=None, aq_params=None):
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

        return reactions, fluxes, bounds, aq_params
    
    @classmethod
    def from_csv(cls, f, bounds=None, aq_params=None):
        raise Exception('ParsedPathway is a virtual class, one should not use its '
                        'static functions.')

    def _to_sbtab(self):
        """
            Returns a full SBtab description of the model.

            Description includes reaction fluxes and per-compound bounds.
            
            This is the superclass implementation, that should be used
            by the inherited classes at the beginning of to_sbtab()
        """
        # Reaction table
        reaction_header = self.SBTAB_GENERIC_HEADER % ('Reaction', 'Reaction')
        reaction_cols = ['!ID', '!ReactionFormula', '!Identifiers:kegg.reaction']

        sio = StringIO()
        sio.writelines([reaction_header + '\n'])
        writer = csv.DictWriter(sio, reaction_cols, dialect='excel-tab')
        writer.writeheader()

        for rxn_id, rxn in zip(self.reaction_ids, self.reactions):
            d = {'!ID': rxn_id,
                 '!ReactionFormula': rxn.GetSlugQueryString(),
                 '!Identifiers:kegg.reaction': rxn.stored_reaction_id}
            writer.writerow(d)

        # Compound table
        compound_header = self.SBTAB_GENERIC_HEADER % ('Compound', 'Compound')
        compound_cols = ['!ID', '!Name', '!Identifiers:kegg.compound']
        sio.writelines(['%\n', compound_header + '\n'])
        writer = csv.DictWriter(sio, compound_cols, dialect='excel-tab')
        writer.writeheader()
        for cid, compound in self.compounds_by_kegg_id.items():
            d = {'!ID': compound.name_slug,
                 '!Name': compound.name.name,
                 '!Identifiers:kegg.compound': cid}
            writer.writerow(d)

        # Flux table
        flux_header = self.SBTAB_GENERIC_HEADER % ('Flux', 'Quantity')
        flux_cols = ['!QuantityType', '!Reaction', '!Reaction:Identifiers:kegg.reaction', '!Value']
        sio.writelines(['%\n', flux_header + '\n'])
        writer = csv.DictWriter(sio, flux_cols, dialect='excel-tab')
        writer.writeheader()

        for rxn_id, rxn, flux in zip(self.reaction_ids, self.reactions, self.fluxes):
            d = {'!QuantityType': 'flux',
                 '!Reaction': rxn_id,
                 '!Reaction:Identifiers:kegg.reaction': rxn.stored_reaction_id,
                 '!Value': flux}
            writer.writerow(d)

        # ConcentrationConstraint table
        conc_header = self.SBTAB_GENERIC_HEADER % ('ConcentrationConstraint', 'Quantity')
        conc_header += " Unit='M'"
        conc_cols = ['!QuantityType', '!Compound',
                     '!Compound:Identifiers:kegg.compound',
                     '!Concentration:Min', '!Concentration:Max']
        sio.writelines(['%\n', conc_header + '\n'])

        writer = csv.DictWriter(sio, conc_cols, dialect='excel-tab')
        writer.writeheader()
        for cid, compound in self.compounds_by_kegg_id.items():
            d = {'!QuantityType': 'concentration',
                 '!Compound': str(compound.name_slug),
                 '!Compound:Identifiers:kegg.compound': cid,
                 '!Concentration:Min': self.bounds.GetLowerBound(cid),
                 '!Concentration:Max': self.bounds.GetUpperBound(cid)}
            writer.writerow(d)

        return sio.getvalue()
    
    def to_sbtab(self):
        raise Exception('ParsedPathway is a virtual class, one should not use its '
                        'static functions.')

    
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

        with sns.axes_style('darkgrid'):
            conc_fig, ax = plt.subplots(1, 1, figsize=(8, 6))
            ax.grid(color='grey', linestyle='--', linewidth=1, alpha=0.2)

            ax.axvspan(1e-8, default_lb, color='y', alpha=0.5)
            ax.axvspan(default_ub, 1e3, color='y', alpha=0.5)
            ax.scatter(concs, ys, label='Variable Concentrations')
            ax.scatter(concs_equal, ys_equal,
                       color='y', label='Fixed Concentrations')
            
            # Special color for metabolites with nonzero shadow prices.
            nz_shadow = np.where(shadow_prices != 0)
            ys_nz_shadow = ys[nz_shadow]
            concs_nz_shadow = concs[nz_shadow]
            ax.scatter(concs_nz_shadow, ys_nz_shadow,
                        color='r', label='Variable Concentrations')
    
            ax.set_yticks(ys)
            ax.set_yticklabels(cnames)
            ax.set_xlabel('Concentration (M)')
            ax.set_xscale('log')
    
            ax.set_xlim(1e-7, 1.0)
            ax.set_ylim(-0.5, Nc - 0.5)
            conc_fig.tight_layout()
    
            svg_data = StringIO()
            conc_fig.savefig(svg_data, format='svg')
        return svg_data.getvalue()
