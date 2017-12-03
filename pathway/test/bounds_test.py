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
        bs = Bounds.from_sbtab(ds[-1])

        for key in bs.lower_bounds:
        	lb = bs.GetLowerBound(key)
        	ub = bs.GetUpperBound(key)
        	msg = 'bounds for %s lb = %.2g, ub = %.2g' % (key, lb, ub)
        	self.assertLessEqual(lb, ub, msg=msg)


if __name__ == '__main__':
    unittest.main()