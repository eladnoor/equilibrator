#!/usr/bin/python

import itertools
import unittest
from util import django_utils

from gibbs import formula_parser


class FormulaParserTest(unittest.TestCase):
    
    def testGetAtomBag(self):
        # (formula, bag)
        test_data = ((None, None),
                     ('C12H22O11', {'C': 12, 'H': 22, 'O': 11}),
                     ('C10H16N5O12P3S', {'C': 10, 'H': 16, 'N': 5,
                                         'O': 12, 'P': 3, 'S': 1}),
                     ('C10R11H16N5O12P3S', {'C': 10, 'H': 16, 'O': 12, 'N': 5, 'P': 3, 'S': 1, 'R': 11}),
                     ('C14H20O4(C5H8)n', {'H': 820, 'C': 514, 'O': 4}),
                     ('C34H32FeN4O4', {'H': 32, 'C': 34, 'Fe': 1, 'O': 4, 'N': 4}))
        
        parser = formula_parser.FormulaParser()
        for formula, expected_bag in test_data:
            #parser.GetAtomBag(formula)
            self.assertEqual(expected_bag, parser.GetAtomBag(formula))
  
  
def Suite():
    suites = (unittest.makeSuite(FormulaParserTest, 'test'))
    return unittest.TestSuite(suites)
    

if __name__ == '__main__':
    unittest.main()