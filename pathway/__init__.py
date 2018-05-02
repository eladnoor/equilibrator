import csv
import logging
import numpy
import seaborn
from io import StringIO

from django.utils.text import slugify
from util import constants
from gibbs import service_config
from gibbs.conditions import AqueousParams
from django.apps import apps
from matplotlib import pyplot as plt
from scipy import linalg
import os
import pandas as pd
from pathway.bounds import Bounds
from pathway.concs import ConcentrationConverter
from pathway.thermo_models import PathwayThermoModel

from equilibrator.settings import BASE_DIR

COFACTORS_FNAME = os.path.join(BASE_DIR, 'pathway/data/cofactors.csv')
DEFAULT_BOUNDS = Bounds.from_csv_filename(
    COFACTORS_FNAME, default_lb=1e-6, default_ub=0.1)


class PathwayParseError(Exception):
    pass


class InvalidReactionFormula(PathwayParseError):
    pass


class UnbalancedReaction(PathwayParseError):
    pass


class ViolatesFirstLaw(PathwayParseError):
    pass


class ParsedPathway(object):
    """A pathway parsed from user input.

    Designed for checking input prior to converting to a stoichiometric model.
    """

    def __init__(self, reactions, fluxes, dG0_r_primes,
                 bounds=None, aq_params=None):
        """Initialize.

        Args:
            reactions: a list of gibbs.reaction.Reaction objects.
            fluxes: numpy.array of relative fluxes in same order as reactions.
            dG0_r_primes: reaction energies.
            bounds: bounds on metabolite concentrations.
                Uses default bounds if None provided.
            aq_params: specify the pH, ionic strength, etc. at which the
                dG values are calculated. May be omitted.
        """
        assert len(reactions) == len(fluxes)
        assert len(reactions) == len(dG0_r_primes)

        self.reactions = reactions
        self.reaction_kegg_ids = [r.stored_reaction_id for r in reactions]
        self.aq_params = aq_params

        self.fluxes = numpy.array(fluxes)
        self.dG0_r_prime = numpy.array(dG0_r_primes)

        self.bounds = bounds or DEFAULT_BOUNDS

        self.S, self.compound_kegg_ids = self._build_stoichiometric_matrix()
        self.compounds_by_kegg_id = self._get_compounds()
        self.compounds = [self.compounds_by_kegg_id[cid]
                          for cid in self.compound_kegg_ids]

        nr, nc = self.S.shape

        # dGr should be orthogonal to nullspace of S
        # If not, dGr is not contained in image(S) and then there
        # is no consistent set of dGfs that generates dGr and the
        # first law of thermo is violated by the model.
        Smat = numpy.matrix(self.S)
        Spinv = linalg.pinv(Smat)
        null_proj = numpy.matrix(numpy.eye(Smat.shape[0])) - Smat*Spinv
        projected = null_proj * numpy.matrix(self.dG0_r_prime).T
        if not numpy.all(projected < 1e-8):
            raise ViolatesFirstLaw(
                'Supplied reaction dG values are inconsistent '
                'with the stoichiometric matrix.')

        # TODO: verify that the vector of standard energies is in the
        # image of the stoichiometric matrix, i.e. that conservation of
        # energy is not violated.

        net_rxn_stoich = (self.fluxes.reshape((nr, 1)) * self.S).sum(axis=0)
        net_rxn_data = []
        for coeff, kid in zip(net_rxn_stoich, self.compound_kegg_ids):
            if coeff != 0:
                net_rxn_data.append(self._reactant_dict(coeff, kid))
        self.net_reaction = apps.get_model('gibbs.reaction').FromIds(net_rxn_data, fetch_db_names=True)
        self._model = self.pathway_model

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
        return ParsedPathway(
            reactions, fluxes, dgs,
            bounds=bounds, aq_params=aq_params)

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
        smat = numpy.zeros((n_reactions, n_compounds))
        for i, s in enumerate(sparses):
            for j, c in enumerate(compounds):
                smat[i, j] = s.get(c, 0)
        return smat, compounds

    @property
    def reactions_balanced(self):
        """Returns true if all pathway reactions are electron and atom-wise balanced."""
        atom_balaned = [r.IsBalanced() for r in self.reactions]
        electron_balaned = [r.IsElectronBalanced() for r in self.reactions]

        balanced = numpy.logical_and(atom_balaned, electron_balaned)
        return numpy.all(balanced)

    @property
    def reactions_have_dG(self):
        return numpy.all([dG is not None for dG in self.dG0_r_prime])

    @property
    def pathway_model(self):
        dGs = numpy.matrix(self.dG0_r_prime).T
        model = PathwayThermoModel(self.S.T, self.fluxes, dGs,
                                   self.compound_kegg_ids,
                                   self.reaction_kegg_ids,
                                   concentration_bounds=self.bounds)
        return model

    def calc_mdf(self):
        model = self.pathway_model
        mdf = model.mdf_result
        return PathwayMDFData(self, mdf)

    def print_reactions(self):
        for f, r in zip(self.fluxes, self.reactions):
            print('%sx %s' % (f, r))

    @classmethod
    def from_full_sbtab(self, reaction_sbtab, flux_sbtab,
                        bounds_sbtab, keqs_sbtab):
        """Returns an initialized ParsedPathway."""
        bounds = Bounds.from_sbtab(bounds_sbtab)

        reaction_df = reaction_sbtab.toDataFrame()
        flux_df = flux_sbtab.toDataFrame()
        bounds_df = bounds_sbtab.toDataFrame()
        keqs_df = keqs_sbtab.toDataFrame()

        name_to_cid = dict(
            zip(bounds_df['Compound'],
                bounds_df['Compound:Identifiers:kegg.compound']))

        query_parser = service_config.Get().query_parser
        reactions = []
        for idx in reaction_df.index:
            row = reaction_df.loc[idx]
            rxn_formula = row['ReactionFormula']
            parsed_rxn = query_parser.ParseReactionQuery(rxn_formula)

            rxn_ds = []
            for coeff, name in parsed_rxn.substrates:
                cid = name_to_cid[name]
                rxn_ds.append(self._reactant_dict(coeff, cid, negate=True))
            for coeff, name in parsed_rxn.products:
                cid = name_to_cid[name]
                rxn_ds.append(self._reactant_dict(coeff, cid, negate=False))
            rxn = apps.get_model('gibbs.reaction').FromIds(rxn_ds, fetch_db_names=True)

            if not rxn.IsBalanced():
                raise UnbalancedReaction(
                    "ReactionFormula '%s' is not balanced" % rxn_formula)
            if not rxn.IsElectronBalanced():
                raise UnbalancedReaction(
                    "ReactionFormula '%s' is not redox balanced" % rxn_formula)

            reactions.append(rxn)

        reaction_ids = reaction_df['ID']
        fluxes = flux_df[flux_df['QuantityType'] == 'flux']
        reaction_fluxes = dict(zip(fluxes['Reaction'], fluxes['Value']))
        fluxes_ordered = [float(reaction_fluxes[rid]) for rid in reaction_ids]

        # grab rows containing keqs.
        keqs = keqs_df[keqs_df['QuantityType'] == 'equilibrium constant']
        reaction_keqs = dict(zip(keqs['Reaction'], keqs['Value']))
        dgs = [-constants.RT * numpy.log(float(reaction_keqs[rid]))
               for rid in reaction_ids]

        # Manually set the delta G values on the reaction objects
        for dg, rxn in zip(dgs, reactions):
            rxn._dg0_prime = dg

        pH = keqs_sbtab.getCustomTableInformation('pH')
        ionic_strength = keqs_sbtab.getCustomTableInformation('IonicStrength')
        ionic_strength_units = keqs_sbtab.getCustomTableInformation(
            'IonicStrengthUnit')
        aq_params = AqueousParams()  # Default values
        if pH:
            aq_params.pH = float(pH)
        if ionic_strength:
            c = float(ionic_strength)
            c = ConcentrationConverter.to_molar_string(
                c, ionic_strength_units)
            aq_params.ionic_strength = c

        pp = ParsedPathway(reactions, fluxes_ordered, dgs,
                           bounds=bounds, aq_params=aq_params)
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
                enz = str(rxn.catalyzing_enzymes[0])
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
            keq = numpy.exp(-dg / constants.RT)
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


class ReactionMDFData(object):

    def __init__(self, reaction, flux, dGr, shadow_price):
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

    @property
    def dG0_prime(self):
        return self.reaction.dg0_prime

    @property
    def dGm_prime(self):
        return self.reaction.dgm_prime


class CompoundMDFData(object):
    def __init__(self, compound, concentration_bounds,
                 concentration, shadow_price):
        self.compound = compound
        self.concentration = concentration
        self.shadow_price = shadow_price
        self.lb, self.ub = concentration_bounds

    @property
    def compound_name(self):
        return self.compound.name.name

    @property
    def link_url(self):
        return '/compound?compoundId=%s' % self.compound.kegg_id

    @property
    def bounds_equal(self):
        return self.lb == self.ub

    @staticmethod
    def html_conc(conc):
        if conc <= 9.999e-4:
            return '%.1f &mu;M' % (1e6*conc)
        return '%.1f mM' % (1e3*conc)

    @property
    def html_concentration(self):
        return self.html_conc(self.concentration)

    @property
    def html_lb(self):
        return self.html_conc(self.lb)

    @property
    def html_ub(self):
        return self.html_conc(self.ub)


class PathwayMDFData(object):

    def __init__(self, parsed_pathway, mdf_result):
        self.parsed_pathway = parsed_pathway
        self.mdf_result = mdf_result
        self.model = mdf_result.model

        rxns = parsed_pathway.reactions
        fluxes = parsed_pathway.fluxes
        dGs = self.mdf_result.dG_r_prime_adj.flatten().tolist()[0]
        prices = self.mdf_result.reaction_prices.flatten().tolist()[0]
        self.reaction_data = [
            ReactionMDFData(*t) for t in zip(rxns, fluxes, dGs, prices)]

        compounds = parsed_pathway.compounds
        cbounds = [self.model.concentration_bounds.GetBoundTuple(cid)
                   for cid in parsed_pathway.compound_kegg_ids]
        concs = self.mdf_result.concentrations.flatten().tolist()[0]
        prices = self.mdf_result.compound_prices.flatten().tolist()[0]
        self.compound_data = [CompoundMDFData(*t) for t in zip(compounds, cbounds, concs, prices)]

    @property
    def mdf(self):
        return self.mdf_result.mdf

    @property
    def min_total_dG(self):
        return self.mdf_result.min_total_dG

    @property
    def max_total_dG(self):
        return self.mdf_result.max_total_dG

    @property
    def max_total_driving_force(self):
        return -self.min_total_dG

    @property
    def min_total_driving_force(self):
        return -self.max_total_dG

    @property
    def conc_plot_svg(self):
        ys = numpy.arange(0, len(self.compound_data))
        concs = numpy.array([c.concentration for c in self.compound_data])
        cnames = [str(c.compound_name) for c in self.compound_data]
        default_lb = self.model.concentration_bounds.default_lb
        default_ub = self.model.concentration_bounds.default_ub

        cids = [str(c.compound.kegg_id) for c in self.compound_data]
        lbs = [self.model.concentration_bounds.GetLowerBound(cid)
               for cid in cids]
        ubs = [self.model.concentration_bounds.GetUpperBound(cid)
               for cid in cids]
        lbs, ubs = numpy.array(lbs), numpy.array(ubs)
        bounds_equal = numpy.where(lbs == ubs)
        ys_equal = ys[bounds_equal]
        concs_equal = concs[bounds_equal]

        # Special color for metabolites with nonzero shadow prices.
        shadow_prices = numpy.array([c.shadow_price for c in self.compound_data])
        nz_shadow = numpy.where(shadow_prices != 0)
        ys_nz_shadow = ys[nz_shadow]
        concs_nz_shadow = concs[nz_shadow]

        conc_figure = plt.figure(figsize=(8, 6))
        seaborn.set_style('darkgrid')
        plt.axes([0.2, 0.1, 0.9, 0.9])
        plt.axvspan(1e-8, default_lb, color='y', alpha=0.5)
        plt.axvspan(default_ub, 1e3, color='y', alpha=0.5)
        plt.scatter(concs, ys, figure=conc_figure,
                    label='Variable Concentrations')
        plt.scatter(concs_equal, ys_equal, figure=conc_figure, color='y',
                    label='Fixed Concentrations')
        plt.scatter(concs_nz_shadow, ys_nz_shadow, figure=conc_figure,
                    color='r', label='Variable Concentrations')

        plt.xticks(family='sans-serif', figure=conc_figure)
        plt.yticks(ys, cnames, family='sans-serif',
            fontsize=8, figure=conc_figure)
        plt.xlabel('Concentration (M)', family='sans-serif',
            figure=conc_figure)
        plt.xscale('log')

        plt.xlim(1e-7, 1.5e2)
        plt.ylim(-1.5, len(self.compound_data) + 0.5)

        svg_data = StringIO()
        conc_figure.savefig(svg_data, format='svg')
        return svg_data.getvalue()

    @property
    def mdf_plot_svg(self):
        dgs = [0] + [r.dGr for r in self.reaction_data]
        dgms = [0] + [r.dGm_prime for r in self.reaction_data]
        cumulative_dgs = numpy.cumsum(dgs)
        cumulative_dgms = numpy.cumsum(dgms)

        xs = numpy.arange(0, len(cumulative_dgs))

        mdf_fig = plt.figure(figsize=(8, 8))
        seaborn.set_style('darkgrid')
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
