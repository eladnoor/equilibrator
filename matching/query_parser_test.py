#!/usr/bin/python
# -*- coding: utf-8 -*-

import query_parser
import unittest


class TestReactionParser(unittest.TestCase):
    """Tests for matcher.Match"""
    
    parsable_reactions = {
        # Simple example
        'H2O => 2 H2 + O2':
        query_parser.ParsedReactionQuery([(1, 'H2O')], [(2, 'H2'), (1, 'O2')]),
        # Two-token names.
        '2 Sodium + 2 Hydrogen chloride <=> 2 Sodium chloride + H2':
        query_parser.ParsedReactionQuery([(2, 'Sodium'), (2, 'Hydrogen chloride')],
                                         [(2, 'Sodium chloride'), (1, 'H2')]),
        # Long reaction, '+' in compound names.
        'Acetyl-CoA + 3 NAD+ + Q + GDP + Pi + 2 H2O <=> '
        'CoA-SH + 3 NADH + 3 H+ + QH2 + GTP + 2 CO2':
        query_parser.ParsedReactionQuery([(1, 'Acetyl-CoA'), (3, 'NAD+'), (1, 'Q'), (1, 'GDP'), (1, 'Pi'), (2, 'H2O')],
                                         [(1, 'CoA-SH'), (3, 'NADH'), (3, 'H+'), (1, 'QH2'), (1, 'GTP'), (2, 'CO2')]),
        # floating point coefficients
        '2 glutamate + 10.5 water -> 6 CO2 + 9 water + NH3': 
        query_parser.ParsedReactionQuery([(2, 'glutamate'), (10.5, 'water')],
                                         [(6, 'CO2'), (9, 'water'), (1, 'NH3')]),
        # fractional coefficients
        'propane + 7/2 o2 => 3 co2 + water':
        query_parser.ParsedReactionQuery([(1, 'propane'), (3.5, 'o2')],
                                         [(3, 'co2'), (1, 'water')]),
        # Names with parens at the beginning
        '(S)-malate => (S)-lactate':
        query_parser.ParsedReactionQuery([(1, '(S)-malate')],
                                         [(1, '(S)-lactate')]),
        'NADH + 1/2 O2 <=> NAD+ + H2O':
        query_parser.ParsedReactionQuery([(1, 'NADH'), (0.5, 'O2')],
                                         [(1, 'NAD+'), (1, 'H2O')]),
        # Unicode arrow thingy
        u'Oxaloacetate + Acetyl-CoA + H2O → Citrate + CoA-SH':
        query_parser.ParsedReactionQuery([(1, 'Oxaloacetate'), (1, 'Acetyl-CoA'), (1, 'H2O')],
                                         [(1, 'Citrate'), (1, 'CoA-SH')])
        }
    
    def setUp(self):
        self._parser = query_parser.QueryParser()
    
    def testIsReactionQuery(self):
        for query in self.parsable_reactions:
            self.assertTrue(self._parser.IsReactionQuery(query),
                            msg='%s is not a reaction query' % query)
    
    def testParseReactions(self):
        for reaction_str, expected_results in self.parsable_reactions.iteritems():
            self.assertTrue(self._parser.IsReactionQuery(reaction_str))
            parsed = self._parser.ParseReactionQuery(reaction_str)
            self.assertEquals(expected_results, parsed)
                
if __name__ == '__main__':
    unittest.main()