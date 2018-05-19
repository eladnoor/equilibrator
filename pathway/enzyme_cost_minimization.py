import sys
from os import path
from io import StringIO
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
from . import ParsedPathway, PathwayAnalysisData, PathwayParseError
import csv
from util import constants
import logging

class EnzymeCostMinimization(ParsedPathway):

    def __init__(self, reactions, fluxes, ecm, bounds=None, aq_params=None):
        super(EnzymeCostMinimization, self).__init__(reactions, fluxes,
             bounds, aq_params)
    
        self.ecm = ecm

    def analyze(self):
        return PathwayECMData(self, self.ecm.ECM())
    
    @classmethod
    def from_sbtab(cls, sbtabs):
        """Returns an initialized ParsedPathway."""
        reactions, fluxes, bounds, aq_params = ParsedPathway._from_sbtab(sbtabs)
        
        # TODO: write a much more detailed model validation function which
        # raises clear exceptions if some data or table is missing!
        
        try:
            sys.path.append(path.expanduser('~/git/enzyme-cost'))
            from ecm import ECMmodel
            ecm = ECMmodel(sbtabs)
        except Exception as e:
            logging.error(str(e))
            raise PathwayParseError('Failed to load the ECM model')
        
        pp = EnzymeCostMinimization(reactions, fluxes, ecm, bounds, aq_params)
        return pp

    @classmethod
    def from_csv(cls, fp, bounds=None, aq_params=None):
        """Returns a pathway parsed from an input file.

        Caller responsible for closing f.

        Args:
            f: file-like object containing CSV data describing the pathway.
        """
        reactions, fluxes, bounds, aq_params = \
            super(EnzymeCostMinimization, cls)._from_csv(fp, bounds, aq_params)
        pp = EnzymeCostMinimization(reactions, fluxes, ecm=None,
                                    bounds=bounds, aq_params=aq_params)
        return pp

    def to_sbtab(self):
        s = self._to_sbtab()
        
        sio = StringIO()
        
        # Parameter table
        # TODO: add all the kinetic parameters (use default values...)
        keq_header = self.SBTAB_GENERIC_HEADER % ('Parameter', 'Quantity')
        keq_cols = ['!QuantityType', '!Reaction', '!Compound', '!Value',
                    '!Unit', '!Reaction:Identifiers:kegg.reaction',
                    '!Compound:Identifiers:kegg.compound', '!ID']
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
        return s + sio.getvalue()

class PathwayECMData(PathwayAnalysisData):
    
    def __init__(self, parsed_pathway, lnC):
        super(PathwayECMData, self).__init__(parsed_pathway)
        self.ecm = parsed_pathway.ecm
        self.lnC = lnC

        enz_concs = self.ecm.ECF(self.lnC)
        for i, r in enumerate(self.reaction_data):
            r.enz_conc = enz_concs[i]

        cid2conc = dict(zip(self.ecm.kegg_model.cids, np.exp(self.lnC)))
        for i, c in enumerate(self.compound_data):
            c.concentration = cid2conc.get(c.compound.kegg_id, 1)
        
    @property
    def score(self):
        """
            Return the total enzyme cost in uM
            (convert from M by multiplying by 10^6)
        """
        return self.ecm.ECF(self.lnC).sum() * 1e6
    
    @property
    def is_ecm(self):
        return True

    @property
    def reaction_plot_svg(self):
        with sns.axes_style('darkgrid'):
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.grid(color='grey', linestyle='--', linewidth=1, alpha=0.2)
            self.ecm.PlotEnzymeDemandBreakdown(self.lnC, ax, plot_measured=False)
            fig.tight_layout()
            
            svg_data = StringIO()
            fig.savefig(svg_data, format='svg')
            return svg_data.getvalue()
