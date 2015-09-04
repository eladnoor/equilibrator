import logging
import matcher
try:
    from nltk.metrics import edit_distance
except ImportError:
    from Levenshtein import distance as edit_distance
import re

from django.db.models import Q
from gibbs import models
from haystack.query import SearchQuerySet


"""
TODO: Delete the Regex and EditDistanceApproxMatchers if the 
haystack-based search works out on the server.
"""


class HaystackApproxMatcher(matcher.Matcher):
    """A matcher that uses the Haystack search plugin."""
    
    def _FindNameMatches(self, query):
        """Override database search."""
        res = SearchQuerySet().autocomplete(name_auto=query)[:2*self._max_results]
        matches = [r.object for r in res]
        return matches
        
    

class RegexApproxMatcher(matcher.Matcher):
    """A matcher that runs regular expressions directly on the database."""
            
    def _PrepareExpression(self, query):
        """Converts the query into a regular expression.
                
        Args:
            query: the string search query.
        
        Returns:
            A regular expression string.
        """
        if not query:
            return None
    
        # Escape regex special characters in the input (search for them
        # literally). Also, we allow '-', ',', '+' and digits in addition to spaces.
        # NOTE(flamholz): We are using MySQL regex syntax here. Might not be 
        # compatible with other databases. See reference:
        #   http://dev.mysql.com/doc/refman/5.1/en/regexp.html
        query = re.escape(query.strip().lower())
        query = re.sub('(\\\?[\s-])+', '[-+,[:digit:][:blank:]]+', query)
        # We allow leading and trailing junk.
        return '.*%s.*' % query

    def _FindNameMatches(self, query):
        """Override database search."""
        expression = self._PrepareExpression(query)
        if not expression:
            return []
        
        matches = models.CommonName.objects.filter(
            name__iregex=expression).prefetch_related(*self._prefetch_objects)
        return matches[:5*self._max_results]
        

class RegexApproxMatcher(matcher.Matcher):
    """A matcher that runs regular expressions directly on the database."""
            
    def _PrepareExpression(self, query):
        """Converts the query into a regular expression.
                
        Args:
            query: the string search query.
        
        Returns:
            A regular expression string.
        """
        if not query:
            return None
    
        # Escape regex special characters in the input (search for them
        # literally). Also, we allow '-', ',', '+' and digits in addition to spaces.
        # NOTE(flamholz): We are using MySQL regex syntax here. Might not be 
        # compatible with other databases. See reference:
        #   http://dev.mysql.com/doc/refman/5.1/en/regexp.html
        query = re.escape(query.strip().lower())
        query = re.sub('(\\\?[\s-])+', '[-+,[:digit:][:blank:]]+', query)
        # We allow leading and trailing junk.
        return '.*%s.*' % query

    def _FindNameMatches(self, query):
        """Override database search."""
        expression = self._PrepareExpression(query)
        if not expression:
            return []
        
        matches = models.CommonName.objects.filter(
            name__iregex=expression).prefetch_related(*self._prefetch_objects)
        return matches[:5*self._max_results]


class EditDistanceMatcher(RegexApproxMatcher):
    """Does near-exact matching and then uses edit distance to score."""
    
    def _GetScore(self, query, match):
        """Custom edit-distance based scoring."""
        str_query = str(query)
        str_candidate = str(match.key)
        dist = float(edit_distance(str_query, str_candidate))
        max_len = float(max(len(str_query), len(str_candidate)))
        return (max_len - dist) / max_len
    
    def _FindNameMatches(self, query):
        """Override database search."""
        qlen = len(query)
        if qlen < 5:
            matches = models.CommonName.objects.filter(
                name__icontains=query).prefetch_related(*self._prefetch_objects)
            return matches[:5*self._max_results]
        
        midpoint = qlen / 2
        head_re = self._PrepareExpression(query[:midpoint])
        tail_re = self._PrepareExpression(query[midpoint:])
        
        matches = models.CommonName.objects.filter(
            Q(name__iregex=head_re) | Q(name__iregex=tail_re)).prefetch_related(*self._prefetch_objects)
        return matches[:5*self._max_results]
    

class CascadingMatcher(matcher.Matcher):
    """A matcher that tries multiple matching strategies."""
    
    def __init__(self, max_results=10, min_score=0.0, match_enzymes=True):
        matcher.Matcher.__init__(self, max_results, min_score, match_enzymes)
        self._exact_matcher = matcher.Matcher(max_results, min_score, match_enzymes)
        self._approx_matcher = HaystackApproxMatcher(max_results, min_score)
    
    def Match(self, query):
        """Override base matching implementation."""  
        matches = self._exact_matcher.Match(query)
        
        match_set = set(m.key for m in matches)
        approx_matches = self._approx_matcher.Match(query)
        for m in approx_matches:
            if m.key not in match_set:
                matches.append(m)
        
        logging.debug('Found %d matches, while the maximum is %d' %
                      (len(matches), self._max_results))
        return self._SortAndClip(matches)