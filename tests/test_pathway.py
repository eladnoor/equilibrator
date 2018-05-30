#!/usr/bin/python3
import os
from django.test import Client
import django
import unittest
import warnings
from equilibrator import settings
import re
import logging
from util.SBtab.SBtabDict import SBtabDict
from pathway import MaxMinDrivingForce, EnzymeCostMinimization, ParsedPathway
from pathway.bounds import Bounds
from pathway.concs import ConcentrationConverter, NoSuchUnits
from gibbs.conditions import AqueousParams
from equilibrator.settings import BASE_DIR
COFACTORS_FNAME = os.path.join(BASE_DIR, 'pathway/data/cofactors.csv')

def ignore_warnings(test_func):
    def do_test(self, *args, **kwargs):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            test_func(self, *args, **kwargs)
    return do_test

class PathwayTester(unittest.TestCase):
    
    def __init__(self, *args, **kwargs):
        super(PathwayTester, self).__init__(*args, **kwargs)
        logging.getLogger().setLevel(logging.WARNING)
        os.environ.setdefault("DJANGO_SETTINGS_MODULE", "equilibrator.settings")
        django.setup()
        self.csv_fname = os.path.join(settings.BASE_DIR, 'tests',
            'example_pathway_ethanol_fermentation.csv')
        self.mdf_sbtab_fname = os.path.join(settings.BASE_DIR, 'tests',
            'example_pathway_ethanol_fermentation_pH7.00_I0.10_MDF.tsv')
        self.ecm_sbtab_fname = os.path.join(settings.BASE_DIR, 'tests',
            'example_pathway_ethanol_fermentation_pH7.00_I0.10_ECM.tsv')

        self.client = Client()

    @ignore_warnings
    def test_mdf_from_csv_file(self):
        aq_params = AqueousParams(pH=7.0, ionic_strength=0.1)
        bounds = Bounds.from_csv_filename(
            COFACTORS_FNAME, default_lb=1e-6, default_ub=1e-2)

        with open(self.csv_fname, 'r') as fp:
            path = MaxMinDrivingForce.from_csv(fp, bounds=bounds,
                                               aq_params=aq_params)
        
        mdf_res = path.analyze()
        self.assertAlmostEqual(mdf_res.score, 1.69, 1)

    @ignore_warnings
    def test_mdf_from_sbtab_file(self):
        with open(self.mdf_sbtab_fname, 'r') as fp:
            sbtabs = SBtabDict.FromSBtabFile(fp)
            
            self.assertSetEqual(set(ParsedPathway.EXPECTED_TNAMES),
                                set(sbtabs.keys()))
            
            bs = Bounds.from_sbtab(sbtabs['ConcentrationConstraint'])
            for key in bs.lower_bounds:
                lb = bs.GetLowerBound(key)
                ub = bs.GetUpperBound(key)
                msg = 'bounds for %s lb = %.2g, ub = %.2g' % (key, lb, ub)
                self.assertLessEqual(lb, ub, msg=msg)

            path = MaxMinDrivingForce.from_sbtab(sbtabs)
            mdf_res = path.analyze()
        
        self.assertAlmostEqual(mdf_res.score, 1.69, 1)
        
    @ignore_warnings
    def test_ecm_from_sbtab_file(self):
        with open(self.ecm_sbtab_fname, 'r') as fp:
            sbtabs = SBtabDict.FromSBtabFile(fp)
            path = EnzymeCostMinimization.from_sbtab(sbtabs)
            ecm_res = path.analyze()
        
        self.assertAlmostEqual(ecm_res.score, 15142.5, 0)

    def test_unit_string(self):
        test_data = [(1.0, 'M', 1.0),
                     (1.0, 'mM', 1e-3),
                     (10, 'mM', 0.01),
                     (100, 'mM', 0.1),
                     (1.0, 'uM', 1e-6),
                     (150, 'millimolar', 0.15),
                     (13.25, 'micromolar', 13.25e-6)]
        for val, unitstring, expected_out in test_data:
            out = ConcentrationConverter.to_molar_string(val, unitstring)
            self.assertEqual(expected_out, out)

    def test_unit_string_failure(self):
        test_data = [(1.0, 'MOLER'),
                     (1.0, 'milimolar'),
                     (1.0, 'kM'),
                     (150, 'kilimoler'),
                     (13.25, 'macromolar')]
        for val, unitstring in test_data:
            self.assertRaises(
                NoSuchUnits,
                ConcentrationConverter.to_molar_string,
                val, unitstring)

    def test_units(self):
        test_data = [(1.0, ConcentrationConverter.UNITS_M, 1.0),
                     (1.0, ConcentrationConverter.UNITS_MM, 1e-3),
                     (1.0, ConcentrationConverter.UNITS_uM, 1e-6),
                     (150, ConcentrationConverter.UNITS_MM, 0.15),
                     (13.25, ConcentrationConverter.UNITS_uM, 13.25e-6)]
        for val, units, expected_out in test_data:
            out = ConcentrationConverter.to_molar_units(val, units)
            self.assertEqual(expected_out, out)

    @ignore_warnings
    def test_web_server(self):
        with open(self.mdf_sbtab_fname, 'r') as fp:
            response = self.client.post('/pathway/results',
                                        {'pathway_file': fp,
                                         'optimization_method': 'MDF'})
        self.assertEqual(response.status_code, 200)
        
        # check that the MDF value is correct
        match = re.findall(r'MDF</a></strong></td>[^<]+<td>([0-9\.]+) kJ/mol</td>',
                           str(response.content))
        for m in match:
            self.assertAlmostEqual(float(m), 1.69, 1)
        
if __name__ == "__main__":
    # PuLP tends to open a TextIOWrapper file without taking care of closing it.
    # Therefore, we have to ignore these warnings, but hopefully in future versions
    # of PuLP that will be solved.
    unittest.main(warnings='ignore')
