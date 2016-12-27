#!/usr/bin/python

import unittest
import logging
from util import django_utils
from util import django_test_utils
from matching import matcher
from django.apps import apps


class TestMatch(unittest.TestCase):
    """Tests for matcher.Match"""

    def testConstruction(self):
        m = matcher.Match('a', True, 0.1)
        self.assertEqual('a', m.key)
        self.assertEqual(True, m.value)
        self.assertEqual(0.1, m.score)


class FirstLastCharacterMatcher(matcher.Matcher):
    """A test matcher.

    Considers a match at the first character perfect and a match
    at the last character slightly worse.
    """

    def _GetScore(self, query, match):
        """Override the single-candidate matching implementation."""
        candidate = str(match.key)
        if query[0] == candidate[0]:
            return 1.0
        if query[-1] == candidate[-1]:
            return 0.75
        return 0.0

    def _FindNameMatches(self, query):
        """Override matching."""
        model = apps.get_model('gibbs.CommonName')
        matches = model.objects.filter(name__icontains=query)
        return matches[:self._max_results]


class TestMatcher(unittest.TestCase):
    """Tests for matcher.Matcher."""

    test_names = ('h2o', 'co2', 'adp', 'nadp', 'abe', 'sol')

    def _CheckIsSortedByScore(self, results):
        prev_score = 10.0  # Scores are all <= 1.0
        for match in results:
            self.assertTrue(match.score <= prev_score)
            prev_score = match.score

    def testBaseMatcher(self):
        logging.info('Testing BaseMatcher')
        m = matcher.Matcher(10)

        for name in self.test_names:
            results = m.Match(name)
            self.assertTrue(len(results) <= 1)

    def testSortingMatcher(self):
        logging.info('Testing SortingMatcher')
        m = FirstLastCharacterMatcher(10)

        results = m.Match('ron')
        self._CheckIsSortedByScore(results)

        results = m.Match('ligand')
        self._CheckIsSortedByScore(results)

    def testMaxResults(self):
        logging.info('Testing MaxResults')
        m = FirstLastCharacterMatcher(2)

        for name in self.test_names:
            results = m.Match(name)
            self.assertTrue(len(results) <= 2)
            self._CheckIsSortedByScore(results)

    def testMinScore(self):
        logging.info('Testing MinScore')
        m = FirstLastCharacterMatcher(max_results=10,
                                      min_score=0.8)

        for name in self.test_names:
            results = m.Match(name)
            self.assertTrue(len(results) <= 10)
            self._CheckIsSortedByScore(results)

            for result in results:
                self.assertTrue(result.score >= 0.8, msg=result.score)


if __name__ == '__main__':
    unittest.main()
