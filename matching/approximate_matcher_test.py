#!/usr/bin/python

from util import django_utils
from util import django_test_utils

import unittest
import approximate_matcher


"""
This test is probably useless now because
1. It doesn't appear to connect to a DB.
2. We are primarily using the haystack based search for matching.
Should figure out how to mock haystack connections for testing.
"""


class TestMatcher(unittest.TestCase):
    """Tests for matcher.Matcher."""

    test_names = ('atp',
                  'adp',
                  'ariel',
                  'glucose',
                  'glucosamine',
                  'alanine',
                  'phenylalanine',
                  'l-glucosamine',
                  'lactate dehydrogenase')

    def _CheckIsSortedByScore(self, results):
        prev_score = 10.0  # Scores are all <= 1.0
        for match in results:
            self.assertTrue(match.score <= prev_score)
            prev_score = match.score

    def _CheckAllNamesOnMatcher(self, names, matcher, max_results, min_score,
                                check_sorted=True):
        for name in names:
            results = matcher.Match(name)
            msg = '%d results but max_results=%d' % (len(results), max_results)
            self.assertTrue(len(results) <= max_results, msg=msg)
            if check_sorted:
                self._CheckIsSortedByScore(results)

            for result in results:
                self.assertTrue(result.score >= min_score)

    def testHaystackMatcher(self, max_results=10, min_score=0.0):
        m = approximate_matcher.HaystackApproxMatcher(
            max_results=max_results, min_score=min_score)
        self._CheckAllNamesOnMatcher(
            self.test_names, m, max_results, min_score, check_sorted=False)

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
