#!/usr/bin/python

import pathways
import unittest

from django.test import TestCase
from util.SBtab import SBtabTools


class PathwayTest(TestCase):

    def test_from_csv_file(self):
    	with open('pathways/test/EMP_glycolysis_simple.csv', 'rU') as f:
    		path = pathways.ParsedPathway.from_csv_file(f)
        	path.print_reactions()

    def test_from_full_sbtab(self):
    	# Happens to be the order in this file.
    	rxns, fluxes, keqs, bounds = SBtabTools.openMultipleSBtab('pathways/test/EMP_glycolysis_full_SBtab.tsv')
    	path = pathways.ParsedPathway.from_full_sbtab(rxns, fluxes, bounds, keqs)
        path.print_reactions()


if __name__ == '__main__':
    unittest.main()