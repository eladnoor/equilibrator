import itertools, logging
from gibbs.models import CommonName, Compound, Enzyme
from haystack.query import SearchQuerySet


class Error(Exception):
    pass


class IllegalQueryError(Error):
    pass


class Match(object):
    """An object containing a string match and it's score."""
    
    def __init__(self, key, value, score):
        """Initialize a Match.
        
        Args:
            key: the object that matched.
            value: the value of the match (the object pointed to by the key).
            score: between 0.0 and 1.0, higher is better.
        """
        self.key = key
        self.value = value
        self.score = score
    
    def __eq__(self, other):
        """Equality checking between matches, used for testing."""
        return (self.key == other.key and
                self.value == other.value and
                self.score == other.score)
    
    def __str__(self):
        """Get as a string for debugging/printability."""
        return '<matcher.Match> value=%s, score=%f' % (self.value,
                                                       self.score)
    
    def TypeStr(self):
        if self.IsCompound():
            return 'Compound'
        elif self.IsEnzyme():
            return 'Enzyme'
        return ''

    def IsCompound(self):
        return isinstance(self.value, Compound)
    
    def IsEnzyme(self):
        return isinstance(self.value, Enzyme)
    
    def Key(self):
        if self.IsCompound():
            return self.value.kegg_id
        elif self.IsEnzyme():
            return self.value.ec
        return None
    

class Matcher(object):
    """A class that matches a string against the database.
    
    The base implementation does exact matching.
    """
    
    def __init__(self, max_results=10, min_score=0.0, match_enzymes=True):
        """Initializes the Matcher.
        
        Args:
            scorer: a MatchScorer object for scoring.
            max_results: the maximum number of matches to return.
            min_score: the minimum match score to return.
        """
        self._max_results = max_results
        self._min_score = min_score
        self._match_enzymes = match_enzymes
        
        self._prefetch_objects = ['compound_set']        
        if self._match_enzymes:
            self._prefetch_objects.extend(['enzyme_set', 'enzyme_set__reactions'])
    
    def _AcceptQuery(self, query):
        """Accept or rejec expression = self._PrepareExpression(query)
        results = models.CommonName.objects.filter(name__iregex=expression)t the query.
        
        Returns:
            True if the query is accepted.
        """
        if query.strip():
            return True
        return False
    
    def _PreprocessQuery(self, query):
        """Perform pre-search query manipulation.
        
        Default implementation simply strips leading/trailing whitespace
        and makes the query lowercase.
        
        Args:
            query: the string query.
            
        Returns:
            The pre-processed query as a string.
        """
        query = query.strip().lower()
        return query
    
    def _PrepocessCandidate(self, candidate):
        """Perform pre-match candidate manipulation.
        
        Default implementation converts to a lower-case string.
        
        Args:
            candidate: the candidate object (convertible to a string).
        
        Returns:
            The pre-processed candidate as a string.
        """
        return str(candidate).strip().lower()
    
    def _FindNameMatches(self, query):
        """Find all the matches for this query.
        
        Args:
            query: the query to match.
            
        Returns:
            A list of CommonName objects matching the query.
        """
        try:
            res = SearchQuerySet().filter(text__exact=query).best_match()
            return [res.object]
        except Exception as e:
            logging.warning('Query failed: ' + str(e))
            return []
    
    def _MakeMatchObjects(self, common_names):
        """Given the list of CommonNames, make the Matches.
        
        Args:
            common_names: a list of CommonNames.
        
        Returns:
            A list of Match objects.
        """
        matches = []
        for name in common_names:
            for compound in name.compound_set.all():
                matches.append(Match(name, compound, 0.0))
            
            if self._match_enzymes:
                for enzyme in name.enzyme_set.all():
                    matches.append(Match(name, enzyme, 0.0))
        
        return matches
        
    def _GetScore(self, query, match):
        """Get the score for a query-match pair.
        
        Args:
            query: the query string.
            match: the Match object.
        
        Returns:
            A score between 0.0 and 1.0.
        """
        query_len = float(len(query))
        candidate_len = float(len(str(match.key)))
        return (query_len / candidate_len)
    
    def _ScoreMatches(self, query, matches):
        """Set the match scores for all matches.
        
        Args:
            query: the query string.
            matches: a list of match objects with uninitialized scores.
        """
        for m in matches:
            m.score = self._GetScore(query, m)
    
    def _FilterMatches(self, matches):
        """Filter the match list for min score.
        
        Args:
            matches: an unfiltered list of match objects.
        """
        # Filter matches without data or beneath the score limit. 
        f = lambda match: (match.score >= self._min_score and
                           match.value)
        filtered = filter(f, matches)
        
        # Take only unique matches.
        group_key = lambda match: match.Key()
        filtered_matches = []
        for _, g in itertools.groupby(filtered, key=group_key):
            # Keep the unique match with the top score.
            max_match = None
            for match in g:
                if not max_match or max_match.score < match.score:
                    max_match = match
            filtered_matches.append(max_match)
            
        return filtered_matches 
    
    def _SortAndClip(self, matches):
        matches.sort(key=lambda m: m.score, reverse=True)
        return matches[:self._max_results]
    
    def Match(self, query):
        """Find matches for the query in the library.
        
        Args:
            query: the string query.
        
        Returns:
            A sorted list of Match objects or None if
            the query could not be parsed.
        """
        if not self._AcceptQuery(query):
            raise IllegalQueryError('%s is not a valid query' % query)
        
        processed_query = self._PreprocessQuery(query)
        logging.debug('Query = %s' % processed_query)
        name_matches = self._FindNameMatches(processed_query)
        logging.debug('Found %d name matches' % len(name_matches))

        matches = self._MakeMatchObjects(name_matches)
        self._ScoreMatches(processed_query, matches)
        matches = self._FilterMatches(matches)
        logging.debug('Found %d matches' % len(matches))
        return self._SortAndClip(matches)
        