#!/usr/bin/python

import unittest
from util import django_utils

# NOTE(flamholz): This is crappy. We're using the real database for
# a unit test. I wish I knew of a better way.
django_utils.SetupDjango()

import approximate_matcher

class TestMatcher(unittest.TestCase):
    """Tests for matcher.Matcher."""
    
    test_names = ('avi',
                  'ari',
                  'ariel',
                  'glucose',
                  'glucosamine',
                  'alanine',
                  'phenylalanine',
                  'l-glucosamine',
                  'lactate dehydrogenase')
    
    def _CheckIsSortedByScore(self, results):
        prev_score = 10.0 # Scores are all <= 1.0
        for match in results:
            self.assertTrue(match.score <= prev_score)
            prev_score = match.score
    
    def _CheckAllNamesOnMatcher(self, names, matcher, max_results, min_score,
                                check_sorted=True):
        for name in names:            
            results = matcher.Match(name)
            self.assertTrue(len(results) <= max_results,
                            msg='%d results but max_results=%d' % (len(results),
                                                                   max_results))
            if check_sorted:
                self._CheckIsSortedByScore(results)
            
            for result in results:
                self.assertTrue(result.score >= min_score)
        
    def testEditDistanceMatcher(self):    
        for max_results in (1, 5, 10):
            for min_score in (0.0, 0.3, 0.7):
                m = approximate_matcher.EditDistanceMatcher(
                    max_results=max_results, min_score=min_score)
                self._CheckAllNamesOnMatcher(
                    self.test_names, m, max_results, min_score)
        

    def testPrepareExpression(self):
        m = approximate_matcher.RegexApproxMatcher({})
        examples = (('  teSt    tEsT ', '.*test[-+,[:digit:][:blank:]]+test.*'),
                    ('gluco', '.*gluco.*'),
                    ('D Fructo', '.*d[-+,[:digit:][:blank:]]+fructo.*'),
                    ('aspartyl-phosphate', '.*aspartyl[-+,[:digit:][:blank:]]+phosphate.*'))
        for query, expression in examples:
            self.assertEqual(expression, m._PrepareExpression(query))

    def testRegexApproxMatcher(self):
        for max_results in (1, 5, 10):
            for min_score in (0.0, 0.3, 0.7):
                m = approximate_matcher.RegexApproxMatcher(
                    max_results=max_results, min_score=min_score)
                self._CheckAllNamesOnMatcher(
                    self.test_names, m, max_results, min_score)
    
    def testCascadingMatcher(self):
        for max_results in (1, 5, 10):
            for min_score in (0.0, 0.3, 0.7):
                m = approximate_matcher.CascadingMatcher(
                    max_results=max_results, min_score=min_score)
                
                # We don't check if things are sorted on the cascading 
                # matcher because it doesn't sort across the results
                # of it's sub-matchers.
                self._CheckAllNamesOnMatcher(
                    self.test_names, m, max_results, min_score,
                    check_sorted=False)


if __name__ == '__main__':
    unittest.main()