#!/usr/bin/python

from pathways.bounds import Bounds
from django.test import TestCase
from util.SBtab import SBtabTools

import logging
import unittest


class TestBounds(TestCase):
    """Tests for matcher.Matcher."""
    
    def testFromSBtab(self):
        ds = SBtabTools.openMultipleSBtab('pathways/test/EMP_glycolysis_full_SBtab.tsv')
        bs = Bounds.from_sbtab(ds[2])
        logging.info(bs.lower_bounds)
        logging.info(bs.upper_bounds)

    
if __name__ == '__main__':
    unittest.main()