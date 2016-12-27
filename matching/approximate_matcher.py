import logging
import matcher
try:
    from nltk.metrics import edit_distance
except ImportError:
    from Levenshtein import distance as edit_distance
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

    def _GetScore(self, query, match):
        """Custom edit-distance based scoring."""
        str_query = str(query)
        str_candidate = str(match.key)
        dist = float(edit_distance(str_query, str_candidate))
        max_len = float(max(len(str_query), len(str_candidate)))
        return (max_len - dist) / max_len

    def _FindNameMatches(self, query):
        """Override database search."""
        # Try plain old autocomplete. If it works, great.
        res = SearchQuerySet().autocomplete(
            title_autocomplete=query)[:self._max_results]
        logging.debug('Found %d results (out of max %d) for "%s" using autocomplete',
                      len(res), self._max_results, query)
        if res:
            return [r.object for r in res]

        # Autocomplete sometimes doesn't work if, for example, you have a
        # spelling error internal to your query. In this case we break the
        # the query into ngrams and search for those ngrams.
        # Sorting is later taken care of by _GetScore.
        res = []
        for i in xrange(len(query) - 3):
            ngram = query[i:i+4]
            auto_res = SearchQuerySet().autocomplete(
                title_autocomplete=ngram)[:self._max_results]
            res.extend(auto_res)
        matches = [r.object for r in res]
        logging.debug('Found %d results using ngrams', len(res))
        return matches


class CascadingMatcher(matcher.Matcher):
    """A matcher that tries multiple matching strategies."""

    def __init__(self, max_results=10, min_score=0.0,
                 match_enzymes=True, return_fast=False):
        matcher.Matcher.__init__(self, max_results, min_score, match_enzymes)
        self._return_fast = return_fast
        self._exact_matcher = matcher.Matcher(
            max_results, min_score, match_enzymes)
        self._approx_matcher = HaystackApproxMatcher(15, min_score)

    def _FilterDups(self, matches):
        """Removes matches pointing to the same primary key."""
        match_ids = set()
        ret = []
        for m in matches:
            mid = '%s:%s', (m.key.TypeStr(), m.value.pk)
            if mid in match_ids:
                # Ignore dups later in the list.
                continue
            match_ids.add(mid)
            ret.append(m)
        return ret

    def Match(self, query):
        """Override base matching implementation."""
        matches = self._exact_matcher.Match(query)
        # In some cases it's advantageous to return exact matches immediately,
        # for example in matching a reaction.
        if matches and self._return_fast:
            logging.debug("Skipping approximate matches for %s", query)
            return self._SortAndClip(matches)

        match_set = set(m.key for m in matches)
        approx_matches = self._approx_matcher.Match(query)
        logging.debug("Approximate matches for %s", query)
        logging.debug([m.key for m in approx_matches])

        for m in approx_matches:
            if m.key not in match_set:
                matches.append(m)

        logging.debug('Found %d matches, while the maximum is %d' %
                      (len(matches), self._max_results))
        matches = self._FilterDups(matches)
        return self._SortAndClip(matches)
