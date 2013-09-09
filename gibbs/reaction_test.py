#!/usr/bin/python

import unittest
from util import django_utils

# NOTE(flamholz): This is crappy. We're using the real database for
# a unit test. I wish I knew of a better way.
django_utils.SetupDjango()

from gibbs import reaction
from gibbs import models


class CompoundWithCoeffTest(unittest.TestCase):
    
    def setUp(self):
        self.compound = models.Compound(kegg_id='fake compound',
                                        formula='C10H16N5O12P3S')
    
    def testMinus(self):
        c_w_coeff = reaction.CompoundWithCoeff(coeff=4, compound=self.compound,
                                               name='foo')
        minus_c = c_w_coeff.Minus()
        self.assertEquals(-4, minus_c.coeff)
        
    def testMicromolarConcentration(self):
        c_w_coeff = reaction.CompoundWithCoeff(coeff=4, compound=self.compound,
                                               name='foo', concentration=0.5)
        self.assertAlmostEqual(5e5, c_w_coeff.micromolar_concentration, 3)


def Suite():
    suites = (unittest.makeSuite(CompoundWithCoeffTest, 'test'),)
    return unittest.TestSuite(suites)

if __name__ == '__main__':
    unittest.main()