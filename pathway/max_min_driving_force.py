import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
import csv
import seaborn as sns
from util import constants
from pathway.thermo_models import PathwayThermoModel
from io import StringIO
from . import ParsedPathway, PathwayAnalyzer, PathwayAnalysisData

class MaxMinDrivingForce(PathwayAnalyzer):
    """
        A class for performing Max-min Driving Force analysis on a given
        pathway (see https://doi.org/10.1371/journal.pcbi.1003483)
    """

    def analyze(self):
        pp = self._parsed_pathway
        model = PathwayThermoModel(pp.S.T, pp.fluxes, pp.dG0_r_primes,
                                   pp.compound_kegg_ids,
                                   pp.reaction_ids,
                                   concentration_bounds=pp.bounds)
        return PathwayMDFData(pp, model.mdf_result)

    @classmethod
    def from_sbtab(cls, sbtabs):
        """
            Reads the input parameters from an SBtabDict object and 
            returns an initialized MaxMinDrivingForce object
        """
        return MaxMinDrivingForce(ParsedPathway.from_sbtab(sbtabs))
    
    @classmethod
    def from_csv(cls, fp, bounds=None, aq_params=None):
        """
            Returns a MaxMinDrivingForce object from an input file.

            Caller responsible for closing fp.
    
            Args:
                f: file-like object containing CSV data describing the pathway.
        """
        parsed_pathway = ParsedPathway.from_csv(fp, bounds, aq_params)
        return MaxMinDrivingForce(parsed_pathway)

    def to_sbtab(self):
        """
            Generates a new SBtab file with all the requirements for running
            MDF. Typically this will be used to convert a CSV file (with
            only the reaction list) into a full model for MDF analysis.
        """
        s = self._parsed_pathway.to_sbtab()
        
        sio = StringIO()
        
        # Parameter table
        keq_header = ParsedPathway.create_sbtab_header('Parameter', 'Quantity',
            pH='%.2f' % self.aq_params.pH, 
            IonicStrength='%.2f' % self.aq_params.ionic_strength,
            IonicStrengthUnit='M')
        sio.write(keq_header)
        keq_cols = ['!QuantityType', '!Reaction', '!Compound', '!Value',
                    '!Unit', '!Reaction:Identifiers:kegg.reaction',
                    '!Compound:Identifiers:kegg.compound', '!ID']

        writer = csv.DictWriter(sio, keq_cols, dialect='excel-tab')
        writer.writeheader()
        for rxn_id, rxn in zip(self.reaction_ids, self.reactions):
            keq_id = 'kEQ_%s' % rxn_id
            keq = np.exp(-rxn.dg0_r_prime / constants.RT)
            d = {'!QuantityType': 'equilibrium constant',
                 '!Reaction': rxn_id,
                 '!Compound': None,
                 '!Value': keq,
                 '!Unit': 'dimensionless',
                 '!Reaction:Identifiers:kegg.reaction': rxn.stored_reaction_id,
                 '!Compound:Identifiers:kegg.compound': None,
                 '!ID': keq_id}
            writer.writerow(d)
        
        sio.write('%\n')
        return s + sio.getvalue()

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

        xticks = np.arange(0, len(cumulative_dgs))-0.5
        xticklabels = [''] + self.reaction_names
        with sns.axes_style('darkgrid'):
            mdf_fig, ax = plt.subplots(1, 1, figsize=(8, 6))
            ax.grid(color=self.COLOR_GRID, linestyle='--',
                    linewidth=1, alpha=0.2)
            ax.plot(cumulative_dgms,
                    label='Characteristic physiological 1 mM concentrations',
                    color=self.COLOR_DELTA_G_M, zorder=1)
            ax.plot(cumulative_dgs,
                    label='MDF-optimized concentrations',
                    color=self.COLOR_DELTA_G_MDF, zorder=1)

            bottleneck_idx = [i for i, r in enumerate(self.reaction_data)
                              if abs(r.shadow_price) != 0]
            lines = [[(i, cumulative_dgs[i]), (i+1, cumulative_dgs[i+1])]
                     for i in bottleneck_idx]
            lines = LineCollection(lines, label='Bottleneck reactions',
                                   linewidth=2,
                                   color=self.COLOR_BOTTLENECK_REACTIONS,
                                   linestyle='-',
                                   zorder=2, alpha=1)
            ax.add_collection(lines)
            ax.set_xticks(xticks)
            ax.set_xticklabels(xticklabels, rotation=45, ha='right')
            ax.set_xlim(0, len(cumulative_dgs)-1)
            mdf_fig.tight_layout()
    
            #ax.set_xlabel('After Reaction Step', family='sans-serif')
            ax.set_xlabel('')
            ax.set_ylabel("Cumulative $\Delta_r G'$ (kJ/mol)", family='sans-serif')
            ax.legend(loc='best', framealpha=0.5)

            svg_data = StringIO()
            mdf_fig.savefig(svg_data, format='svg')
        
        return svg_data.getvalue()

