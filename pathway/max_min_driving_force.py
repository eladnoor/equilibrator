import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
import csv
import seaborn as sns
from util import constants
from pathway.thermo_models import PathwayThermoModel
from io import StringIO
from . import ParsedPathway, PathwayAnalysisData

class MaxMinDrivingForce(ParsedPathway):
    """A pathway parsed from user input.

    Designed for checking input prior to converting to a stoichiometric model.
    """

    def analyze(self):
        model = PathwayThermoModel(self.S.T, self.fluxes, self.dG0_r_primes,
                                   self.compound_kegg_ids,
                                   self.reaction_ids,
                                   concentration_bounds=self.bounds)
        return PathwayMDFData(self, model.mdf_result)

    @classmethod
    def from_sbtab(cls, sbtabs):
        """Returns an initialized ParsedPathway."""
        pp = MaxMinDrivingForce(*ParsedPathway._from_sbtab(sbtabs))
        return pp
    
    @classmethod
    def from_csv(cls, fp, bounds=None, aq_params=None):
        """Returns a pathway parsed from an input file.

        Caller responsible for closing f.

        Args:
            f: file-like object containing CSV data describing the pathway.
        """
        reactions, fluxes, bounds, aq_params = \
            super(MaxMinDrivingForce, cls)._from_csv(fp, bounds, aq_params)
        pp = MaxMinDrivingForce(reactions, fluxes,
                                bounds=bounds, aq_params=aq_params)
        return pp

    def to_sbtab(self):
        s = self._to_sbtab()
        
        sio = StringIO()
        
        # Parameter table
        keq_header = self.create_sbtab_header('Parameter', 'Quantity',
            pH='%.2f' % self.aq_params.pH, 
            IonicStrength='%.2f' % self.aq_params.ionic_strength,
            IonicStrengthUnit='M')

        keq_cols = ['!QuantityType', '!Reaction', '!Compound', '!Value',
                    '!Unit', '!Reaction:Identifiers:kegg.reaction',
                    '!Compound:Identifiers:kegg.compound', '!ID']
        sio.writelines(['%\n', keq_header + '\n'])

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
            ax.grid(color=ParsedPathway.COLOR_GRID, linestyle='--',
                    linewidth=1, alpha=0.2)
            ax.plot(cumulative_dgms,
                    label='Characteristic physiological 1 mM concentrations',
                    color=ParsedPathway.COLOR_DELTA_G_M, zorder=1)
            ax.plot(cumulative_dgs,
                    label='MDF-optimized concentrations',
                    color=ParsedPathway.COLOR_DELTA_G_MDF, zorder=1)

            bottleneck_idx = [i for i, r in enumerate(self.reaction_data)
                              if abs(r.shadow_price) != 0]
            lines = [[(i, cumulative_dgs[i]), (i+1, cumulative_dgs[i+1])]
                     for i in bottleneck_idx]
            lines = LineCollection(lines, label='Bottleneck reactions',
                                   linewidth=2,
                                   color=ParsedPathway.COLOR_BOTTLENECK_REACTIONS,
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

