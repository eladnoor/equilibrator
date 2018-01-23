import logging
from matching import matcher
from haystack.query import SearchQuerySet


class HaystackApproxMatcher(matcher.Matcher):
    """A matcher that uses the Haystack search plugin.

    Current behavior:
        First uses the haystack autocomplete. If there are results, returns.
        If no results, then there are no exact matches for your search,
        in which case we break the query into 4-grams and search for those.
        We then let the parent class logic dedup those and they are ranked
        according to their edit-distance to the query, as per _GetScore below.
    """

    def _FindNameMatches(self, query):
        """Override database search."""
        # Try plain old autocomplete. If it works, great.
        logging.debug('Trying Autocomplete for query "%s"', query)
        res = SearchQuerySet().autocomplete(title_autocomplete=query)
        
        if len(res) > self._max_results:
            logging.debug('Autocomplete found %d matches for "%s", keeping first %d',
                          len(res), query, self._max_results)
            return [r.object for r in res[0:self._max_results]]
        elif len(res) > 0:
            logging.debug('Autocomplete found %d matches for "%s"',
                          len(res), query)
            return [r.object for r in res]
        
        logging.debug('Autocomplete found no matches for "%s", trying N-grams',
                      query)
        # Autocomplete sometimes doesn't work if, for example, you have a
        # spelling error internal to your query. In this case we break the
        # the query into ngrams and search for those ngrams.
        # Sorting is later taken care of by _GetScore.
        res = []
        for i in range(len(query) - 3):
            ngram = query[i:i+4]
            auto_res = SearchQuerySet().autocomplete(
                title_autocomplete=ngram)[:self._max_results]
            res.extend(auto_res)
        logging.debug('N-grams found %d matches for "%s"', len(res), query)
        return [r.object for r in res]

class CascadingMatcher(matcher.Matcher):
    """A matcher that tries multiple matching strategies."""

    def __init__(self, max_results=10, min_score=0.0,
                 match_enzymes=True, return_fast=False):
        matcher.Matcher.__init__(self, max_results, min_score, match_enzymes)
        self._return_fast = return_fast
        self._exact_matcher = matcher.Matcher(
            max_results, min_score, match_enzymes)
        self._approx_matcher = HaystackApproxMatcher(15, min_score)

    def Match(self, query):
        """Override base matching implementation."""
        matches = self._exact_matcher.Match(query)

        # In some cases it's advantageous to return exact matches immediately,
        # for example in matching a reaction.
        if matches and self._return_fast:
            logging.debug("Skipping approximate matches for %s", query)
            return self._SortAndClip(matches)

        logging.debug("Approximate matches for %s", query)
        approx_matches = self._approx_matcher.Match(query)
        matches += approx_matches
        matches = self._FilterMatches(matches)
        matches = self._SortAndClip(matches)
        return matches
