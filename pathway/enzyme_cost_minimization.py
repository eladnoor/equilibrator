import sys
from os import path
from io import StringIO
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
from . import ParsedPathway, PathwayAnalysisData

class EnzymeCostMinimization(ParsedPathway):

    EXPECTED_TNAMES = ['Reaction', 'Compound', 'Parameter', 'Flux',
                       'ConcentrationConstraint']

    def __init__(self, reactions, fluxes, bounds=None, aq_params=None):
        super(EnzymeCostMinimization, self).__init__(reactions, fluxes,
             bounds, aq_params)
    
        self.ecm = None

    def is_empty(self):
        return self.ecm.ecf.Nr == 0
    
    def analyze(self):
        lnC = self.ecm.ECM()
        return PathwayECMData(self, lnC)
    
    @classmethod
    def from_sbtabs(cls, sbtabs):
        """Returns an initialized ParsedPathway."""
        sys.path.append(path.expanduser('~/git/enzyme-cost'))
        from ecm import ECMmodel

        reactions, fluxes, bounds = \
            ParsedPathway.get_data_from_sbtabs(sbtabs)
        
        pp = EnzymeCostMinimization(reactions, fluxes, bounds)
        pp.ecm = ECMmodel(sbtabs)
        return pp


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
        fig, ax = plt.subplots(figsize=(8, 8))
        sns.set_style('darkgrid')
        self.ecm.PlotEnzymeDemandBreakdown(self.lnC, ax, plot_measured=False)

        svg_data = StringIO()
        fig.savefig(svg_data, format='svg')
        return svg_data.getvalue()
