#!/usr/bin/python

import unittest
from util.SBtab import SBtabTools


class SBtabTest(unittest.TestCase):
    
    def testOpenMultipleSBtab(self):
        ds = SBtabTools.openMultipleSBtab('pathways/test/EMP_glycolysis_full_SBtab.tsv')
        self.assertEqual(3, len(ds))

    def testOpenMultipleSBtabFromFile(self):
		with open('pathways/test/EMP_glycolysis_full_SBtab.tsv') as f:
			ds = SBtabTools.openMultipleSBtabFromFile(f)
			self.assertEqual(3, len(ds))


if __name__ == '__main__':
    unittest.main()
