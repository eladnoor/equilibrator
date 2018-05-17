import numpy as np
from matplotlib import pyplot as plt
import csv
import logging
import seaborn as sns
from django.apps import apps
from util import constants
from gibbs import service_config
from gibbs.conditions import AqueousParams
import pandas as pd
from pathway.thermo_models import PathwayThermoModel
from django.utils.text import slugify
from pathway.concs import ConcentrationConverter
from io import StringIO
from scipy import linalg
from . import InvalidReactionFormula, UnbalancedReaction, ViolatesFirstLaw, \
    ParsedPathway, PathwayAnalysisData, ReactionData, CompoundData

class MaxMinDrivingForce(ParsedPathway):
    """A pathway parsed from user input.

    Designed for checking input prior to converting to a stoichiometric model.
    """

    EXPECTED_TNAMES = ['Reaction', 'RelativeFlux', 'ReactionConstant',
                       'ConcentrationConstraint']

    def __init__(self, reactions, fluxes, bounds=None, aq_params=None):
        super(MaxMinDrivingForce, self).__init__(reactions, fluxes,
             bounds, aq_params)
    
        self.dG0_r_prime = None

    @classmethod
    def from_csv_file(cls, f,
                      bounds=None, aq_params=None):
        """Returns a pathway parsed from an input file.

        Caller responsible for closing f.

        Args:
            f: file-like object containing CSV data describing the pathway.
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

            reactions.append(rxn)

        dgs = [r.DeltaG0Prime(aq_params) for r in reactions]
        pp = MaxMinDrivingForce(reactions, fluxes,
                                bounds=bounds, aq_params=aq_params)
        pp.dG0_r_prime = dgs
        return pp

    @property
    def reactions_balanced(self):
        """Returns true if all pathway reactions are electron and atom-wise balanced."""
        atom_balaned = [r.IsBalanced() for r in self.reactions]
        electron_balaned = [r.IsElectronBalanced() for r in self.reactions]

        balanced = np.logical_and(atom_balaned, electron_balaned)
        return np.all(balanced)

    @property
    def reactions_have_dG(self):
        return np.all([dG is not None for dG in self.dG0_r_prime])

    @property
    def pathway_model(self):
        dGs = np.matrix(self.dG0_r_prime).T
        model = PathwayThermoModel(self.S.T, self.fluxes, dGs,
                                   self.compound_kegg_ids,
                                   self.reaction_kegg_ids,
                                   concentration_bounds=self.bounds)
        return model

    def analyze(self):
        return PathwayMDFData(self, self.pathway_model.mdf_result)

    def print_reactions(self):
        for f, r in zip(self.fluxes, self.reactions):
            print('%sx %s' % (f, r))

    @classmethod
    def from_sbtabs(cls, sbtabs):
        """Returns an initialized ParsedPathway."""
        reactions, fluxes, bounds = \
            ParsedPathway.get_data_from_sbtabs(sbtabs)
        
        reaction_ids = sbtabs['Reaction'].toDataFrame()['ID']

        # read the general aqueous parameters from the ReactionConstant
        # table header
        reaction_constant_sbtab = sbtabs['ReactionConstant']
        pH, ionic_strength, ionic_strength_units = \
            map(reaction_constant_sbtab.getCustomTableInformation,
                ['pH', 'IonicStrength', 'IonicStrengthUnit'])
        aq_params = AqueousParams()  # Default values
        if pH:
            aq_params.pH = float(pH)
        if ionic_strength:
            c = float(ionic_strength)
            c = ConcentrationConverter.to_molar_string(c, ionic_strength_units)
            aq_params.ionic_strength = c

        # create the ParsedPathway object
        pp = MaxMinDrivingForce(reactions, fluxes, bounds, aq_params)
        
        # read Keq table
        keqs_df = reaction_constant_sbtab.toDataFrame()
        keqs = keqs_df[keqs_df['QuantityType'] == 'equilibrium constant']
        reaction_keqs = dict(zip(keqs['Reaction'], keqs['Value']))
        dgs = [-constants.RT * np.log(float(reaction_keqs[rid]))
               for rid in reaction_ids]

        # Manually set the delta G values on the reaction objects
        for dg, rxn in zip(dgs, pp.reactions):
            rxn._dg0_prime = dg

        # dGr should be orthogonal to nullspace of S
        # If not, dGr is not contained in image(S) and then there
        # is no consistent set of dGfs that generates dGr and the
        # first law of thermo is violated by the model.
        pp.dG0_r_prime = np.array(dgs)
        Smat = np.matrix(pp.S)
        Spinv = linalg.pinv(Smat)
        null_proj = np.matrix(np.eye(Smat.shape[0])) - Smat*Spinv
        projected = null_proj * np.matrix(pp.dG0_r_prime).T
        if not np.all(projected < 1e-8):
            raise ViolatesFirstLaw(
                'Supplied reaction dG values are inconsistent '
                'with the stoichiometric matrix.')

        return pp

    def to_full_sbtab(self):
        """Returns a full SBtab description of the model.

        Description includes reaction fluxes and per-compound bounds.
        """
        generic_header_fmt = "!!SBtab TableName='%s' TableType='%s' Document='%s' SBtabVersion='1.0'"
        reaction_header = generic_header_fmt % ('Reaction', 'Reaction', 'Pathway Model')
        reaction_cols = ['!ID', '!ReactionFormula', '!Identifiers:kegg.reaction']

        sio = StringIO()
        sio.writelines([reaction_header + '\n'])
        writer = csv.DictWriter(sio, reaction_cols, dialect='excel-tab')
        writer.writeheader()

        rxn_ids = []
        for i, rxn in enumerate(self.reactions):
            kegg_id = rxn.stored_reaction_id
            rxn_id = kegg_id
            if rxn.catalyzing_enzymes:
                enz = str(rxn.catalyzing_enzymes[0].FirstName().name)
                enz_slug = slugify(enz)[:10]
                enz_slug = enz_slug.replace('-', '_')
                rxn_id = '%s_%s' % (enz_slug, kegg_id)
            elif not kegg_id:
                rxn_id = 'RXN%03d' % i

            rxn_ids.append(rxn_id)
            d = {'!ID': rxn_id,
                 '!ReactionFormula': rxn.GetSlugQueryString(),
                 '!Identifiers:kegg.reaction': kegg_id}
            writer.writerow(d)

        # Relative fluxes
        flux_header = generic_header_fmt % ('RelativeFlux', 'Quantity', 'Pathway Model')
        flux_cols = ['!QuantityType', '!Reaction', '!Reaction:Identifiers:kegg.reaction', '!Value']
        sio.writelines(['%\n', flux_header + '\n'])
        writer = csv.DictWriter(sio, flux_cols, dialect='excel-tab')
        writer.writeheader()

        for i, rxn_id in enumerate(rxn_ids):
            d = {'!QuantityType': 'flux',
                 '!Reaction': rxn_id,
                 '!Reaction:Identifiers:kegg.reaction': self.reactions[i].stored_reaction_id,
                 '!Value': self.fluxes[i]}
            writer.writerow(d)

        # Write KEQs.
        keq_header = generic_header_fmt % (
            'ReactionConstant', 'Quantity', 'Pathway Model')
        keq_cols = ['!QuantityType', '!Reaction', '!Value',
                    '!Unit', '!Reaction:Identifiers:kegg.reaction', '!ID']
        if self.aq_params:
            # Write pH and ionic strength in header
            aq_params_header = (
                "pH='%.2f' 'IonicStrength='%.2f' IonicStrengthUnit='M'")

            aq_params_header = aq_params_header % (
                self.aq_params.pH, self.aq_params.ionic_strength)
            keq_header = '%s %s' % (keq_header, aq_params_header)
        sio.writelines(['%\n', keq_header + '\n'])

        writer = csv.DictWriter(sio, keq_cols, dialect='excel-tab')
        writer.writeheader()
        for i, (rxn_id, rxn, dg) in enumerate(zip(rxn_ids, self.reactions, self.dG0_r_prime)):
            keq_id = 'kEQ_R%d' % i
            keq = np.exp(-dg / constants.RT)
            d = {'!QuantityType': 'equilibrium constant',
                 '!Reaction': rxn_id,
                 '!Value': keq,
                 '!Unit': 'dimensionless',
                 '!Reaction:Identifiers:kegg.reaction': rxn.stored_reaction_id,
                 '!ID': keq_id}
            writer.writerow(d)

        conc_header = generic_header_fmt % ('ConcentrationConstraint', 'Quantity', 'Pathway Model')
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


class PathwayMDFData(PathwayAnalysisData):

    def __init__(self, parsed_pathway, mdf_result):
        super(PathwayMDFData, self).__init__(parsed_pathway)
        self.mdf_result = mdf_result

        for i, r in enumerate(self.reaction_data):
            r.dGr = self.mdf_result.dG_r_prime_adj[i, 0]
            r.shadow_price = self.mdf_result.reaction_prices[i, 0]

        for i, c in enumerate(self.compound_data):
            c.concentration = self.mdf_result.concentrations[i, 0]
            c.shadow_price = self.mdf_result.compound_prices[i, 0]

    @property
    def score(self):
        return self.mdf_result.mdf

    @property
    def is_mdf(self):
        return True

    @property
    def max_total_driving_force(self):
        return -self.mdf_result.min_total_dG

    @property
    def min_total_driving_force(self):
        return -self.mdf_result.max_total_dG

    @property
    def reaction_plot_svg(self):
        dgs = [0] + [r.dGr for r in self.reaction_data]
        dgms = [0] + [r.dGm_prime for r in self.reaction_data]
        cumulative_dgs = np.cumsum(dgs)
        cumulative_dgms = np.cumsum(dgms)

        xs = np.arange(0, len(cumulative_dgs))

        mdf_fig = plt.figure(figsize=(8, 8))
        sns.set_style('darkgrid')
        plt.plot(xs, cumulative_dgms,
                 label='Characteristic physiological 1 mM concentrations')
        plt.plot(xs, cumulative_dgs,
                 label='MDF-optimized concentrations')
        plt.xticks(xs, family='sans-serif')
        plt.yticks(family='sans-serif')

        # TODO: Consider using reaction IDs from the file as xticks?
        plt.xlabel('After Reaction Step', family='sans-serif')
        plt.ylabel("Cumulative $\Delta_r G'$ (kJ/mol)", family='sans-serif')
        plt.legend(loc=3)

        svg_data = StringIO()
        mdf_fig.savefig(svg_data, format='svg')
        return svg_data.getvalue()

