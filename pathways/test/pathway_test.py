#!/usr/bin/python

import pathways
import unittest

from django.test import TestCase
from util.SBtab import SBtabTools


class PathwayTest(TestCase):

    def test_from_csv_file(self):
    	with open('pathways/test/EMP_glycolysis_to_EtOH_simple.tsv') as f:
    		path = pathways.ParsedPathway.from_csv_file(f)
        	path.print_reactions()

    def _notest_from_simple_sbtab(self):
    	ds = SBtabTools.openMultipleSBtab('pathways/test/EMP_glycolysis_to_EtOH_simple.tsv')
    	path = pathways.ParsedPathway.from_simple_sbtab(ds[0])
        path.print_reactions()


if __name__ == '__main__':
    unittest.main()