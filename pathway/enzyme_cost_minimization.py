import sys
from os import path
from io import StringIO
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
from . import ParsedPathway, PathwayAnalyzer, PathwayAnalysisData, PathwayParseError
import csv
from util import constants
import logging

try:
    sys.path.append(path.expanduser('~/git/enzyme-cost'))
    from ecm import ECMmodel
except ImportError as e:
    logging.error('Make sure to install the enzyme-cost package')
    raise e

class EnzymeCostMinimization(PathwayAnalyzer):
    """
        A class for performing Enzyme Cost Minimization analysis on a given
        pathway (see https://doi.org/10.1371/journal.pcbi.1005167)
    """

    def __init__(self, parsed_pathway, ecm):
        super(EnzymeCostMinimization, self).__init__(parsed_pathway)
        self.ecm = ecm

    def analyze(self):
        return PathwayECMData(self._parsed_pathway, self.ecm, self.ecm.ECM())
    
    @classmethod
    def from_sbtab(cls, sbtabs):
        """Returns an initialized ParsedPathway."""
        parsed_pathway = ParsedPathway.from_sbtab(sbtabs)
        
        try:
            ecm = ECMmodel(sbtabs)
        except Exception as e:
            logging.error(str(e))
            raise PathwayParseError('Failed to load the ECM model')
            
        # TODO: write a much more detailed model validation function which
        # raises clear exceptions if some data or table is missing!
        
        return EnzymeCostMinimization(parsed_pathway, ecm)

    @classmethod
    def from_csv(cls, fp, bounds=None, aq_params=None):
        """Returns an EnzymeCostMinimization from an input file.

        Caller responsible for closing f.

        Args:
            f: file-like object containing CSV data describing the pathway.
        """
        parsed_pathway = ParsedPathway.from_csv(fp, bounds, aq_params)
        ecm = None
        return EnzymeCostMinimization(parsed_pathway, ecm)

    def to_sbtab(self):
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

        # equilibrium constants        
        for rxn_id, rxn in zip(self.reaction_ids, self.reactions):
            param_id = 'kEQ_%s' % rxn_id
            keq = np.exp(-rxn.dg0_r_prime / constants.RT)
            d = {'!QuantityType': 'equilibrium constant',
                 '!Reaction': rxn_id,
                 '!Compound': None,
                 '!Value': keq,
                 '!Unit': 'dimensionless',
                 '!Reaction:Identifiers:kegg.reaction': rxn.stored_reaction_id,
                 '!Compound:Identifiers:kegg.compound': None,
                 '!ID': param_id}
            writer.writerow(d)
            
        # add default values for the 3 types catalytic rate constants
        qtype_list = [('catalytic rate constant geometric mean', 'kC'),
                      ('substrate catalytic rate constant', 'kcrf'),
                      ('product catalytic rate constant', 'kcrr')]
        
        for qtype, id_prefix in qtype_list:
            for rxn_id, rxn in zip(self.reaction_ids, self.reactions):
                param_id = '%s_%s' % (id_prefix, rxn_id)
                d = {'!QuantityType': qtype,
                     '!Reaction': rxn_id,
                     '!Compound': None,
                     '!Value': 100.0,
                     '!Unit': '1/s',
                     '!Reaction:Identifiers:kegg.reaction': rxn.stored_reaction_id,
                     '!Compound:Identifiers:kegg.compound': None,
                     '!ID': param_id}
                writer.writerow(d)

        for rxn_id, rxn in zip(self.reaction_ids, self.reactions):
            for cpd in map(lambda c: c.compound, rxn.substrates + rxn.products):
                param_id = 'kM_%s_%s' % (rxn_id, cpd.name_slug)
                d = {'!QuantityType': 'Michaelis constant',
                     '!Reaction': rxn_id,
                     '!Compound': cpd.name_slug,
                     '!Value': 1.0,
                     '!Unit': 'mM',
                     '!Reaction:Identifiers:kegg.reaction': rxn.stored_reaction_id,
                     '!Compound:Identifiers:kegg.compound': cpd.kegg_id,
                     '!ID': param_id}
                writer.writerow(d)                
            
        # protein molecular mass            
        for rxn_id, rxn in zip(self.reaction_ids, self.reactions):
            d = {'!QuantityType': 'protein molecular mass',
                 '!Reaction': rxn_id,
                 '!Compound': None,
                 '!Value': 10000,
                 '!Unit': 'Da',
                 '!Reaction:Identifiers:kegg.reaction': rxn.stored_reaction_id,
                 '!Compound:Identifiers:kegg.compound': None,
                 '!ID': None}
            writer.writerow(d)
        
        # compound moleulcar mass
        for cpd in self.compounds:
            d = {'!QuantityType': 'molecular mass',
                 '!Reaction': None,
                 '!Compound': cpd.name_slug,
                 '!Value': '%.1f' % cpd.mass,
                 '!Unit': 'Da',
                 '!Reaction:Identifiers:kegg.reaction': None,
                 '!Compound:Identifiers:kegg.compound': cpd.kegg_id,
                 '!ID': None}
            writer.writerow(d)        

        sio.write('%\n')

        return s + sio.getvalue()

class PathwayECMData(PathwayAnalysisData):
    
    def __init__(self, parsed_pathway, ecm, lnC):
        super(PathwayECMData, self).__init__(parsed_pathway)
        self.ecm = ecm
        self.lnC = lnC

        for rxn, c in zip(self.reaction_data, self.enzyme_costs):
            rxn.enz_conc = c
        
        enz_cost_brkdwn = self.enzyme_costs_breakdown
        for i, rxn in enumerate(self.reaction_data):
            rxn.min_enz_conc = enz_cost_brkdwn[i, 0]
            rxn.eta_th = 1.0/enz_cost_brkdwn[i, 1]
            rxn.eta_sat = 1.0/enz_cost_brkdwn[i, 2]

        cid2conc = dict(zip(self.ecm.kegg_model.cids, np.exp(self.lnC)))
        for i, cpd in enumerate(self.compound_data):
            cpd.concentration = cid2conc.get(cpd.compound.kegg_id, 1)
        
    @property
    def score(self):
        """
            Return the total enzyme cost in uM
            (convert from M by multiplying by 10^6)
        """
        return self.enzyme_costs.sum() * 1e6
    
    @property
    def enzyme_costs(self):
        return self.ecm.ECF(self.lnC)

    @property
    def enzyme_costs_breakdown(self):
        return np.array(self.ecm.ecf.GetEnzymeCostPartitions(self.lnC))
    
    @property
    def is_ecm(self):
        return True

    @property
    def reaction_plot_svg(self):
        """
            Produces a horizontal stacked bar plot of the enzyme demands, broken down
            into 3 parts:
                - ideal demand (i.e. the minimum amount derived from the flux and kcat)
                - inverse of thermodynamic efficiency (1 / eta^th)
                - inverse of saturation efficiency (1 / eta^sat)
            The values are shown in log-scale, since the total demand is the product
            of these 3 values.
        """
        with sns.axes_style('darkgrid'):
            fig, ax = plt.subplots(figsize=(8, 8))
            ax.grid(color='grey', linestyle='--', linewidth=1, alpha=0.2)
            
            labels = self.ecm.ecf.ECF_LEVEL_NAMES[0:3]
            costs = self.enzyme_costs_breakdown[::-1, :] # reverse order of reactions
            
            base = min(filter(None, costs[:, 0])) / 2.0
            idx_zero = (costs[:, 0] == 0)
            costs[idx_zero, 0] = base
            costs[idx_zero, 1:] = 1.0
    
            lefts = np.hstack([np.ones((costs.shape[0], 1)) * base,
                               np.cumprod(costs, 1)])
            steps = np.diff(lefts)
    
            ind = range(costs.shape[0])    # the x locations for the groups
            height = 0.8
            ax.set_xscale('log')
    
            colors = sns.color_palette('muted', n_colors=6)[3:6]
    
            for i, label in enumerate(labels):
                ax.barh(ind, width=steps[:, i].flat, height=height,
                        left=lefts[:, i].flat, color=colors[i])
    
            ax.set_yticks(ind)
            ax.set_yticklabels(reversed(self.reaction_names), size='medium')
            ax.legend(labels, loc='best', framealpha=0.5)
            ax.set_xlabel('enzyme demand [M]')
            ax.set_xlim(base, None)
            fig.tight_layout()
            
            svg_data = StringIO()
            fig.savefig(svg_data, format='svg')
            return svg_data.getvalue()

