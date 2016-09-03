#!/usr/bin/python

from pathways.concs import ConcentrationConverter, NoSuchUnits
from django.test import TestCase
from util.SBtab import SBtabTools

import logging
import unittest


class TestConcentrationConverter(TestCase):
    """Tests for ConcentrationConverter."""
   
    def test_unit_string(self):
        test_data = [(1.0, 'M', 1.0),
                     (1.0, 'mM', 1e-3),
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


if __name__ == '__main__':
    unittest.main()