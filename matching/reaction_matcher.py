import logging

from gibbs import models

class ReactionCompoundMatch(object):
    """A match of a compound in a reaction.
    
    Contains the parsed name (what the user entered), the parsed
    stoichiometric coefficient, and a list of matcher.Match objects
    of the potential matches.
    """
    def __init__(self, parsed_name, parsed_coeff, matches):
        self.parsed_name = parsed_name
        self.parsed_coeff = parsed_coeff
        self.matches = matches
    
    def __str__(self):
        return '%d %s, matches: %s' % (self.parsed_coeff,
                                       self.parsed_name,
                                       ', '.join(self.matches))
    
    def ParsedDataEqual(self, other):
        """Checks if the parsed data (name, coefficient) are equal
           for two ReactionCompoundMathes.
        """
        return (self.parsed_coeff == other.parsed_coeff and
                self.parsed_name == other.parsed_name)

    
class ReactionMatches(object):
    """A reaction parsed from a query with all possible matches."""
    
    def __init__(self, substrates=None, products=None):
        """Initialize the ParsedReaction object.
        
        Args:
            reactants: a list of ReactionCompoundMatches for the reactants.
            products: a list of ReactionCompoundMatches for the products.
        """
        self.substrates = substrates or []
        self.products = products or []
    
    @staticmethod
    def _FindFirstCompoundMatch(matches):
        """Finds the first match that has a Compound as a value."""
        for m in matches:
            if isinstance(m.value, models.Compound):
                return m
        return None
    
    def GetBestMatch(self):
        """Returns a 2-tuple of product and reactant lists for the
        best match.
        
        Each list is of 3 tuples (coeff, kegg_id, name).
        """
        reactants = []
        for c in self.substrates:
            compound_match = self._FindFirstCompoundMatch(c.matches)
            if not compound_match:
                return None
            
            reactants.append((c.parsed_coeff,
                              compound_match.value.kegg_id,
                              compound_match.key))
        
        products = []
        for c in self.products:
            compound_match = self._FindFirstCompoundMatch(c.matches)
            if not compound_match:
                return None
            
            products.append((c.parsed_coeff,
                             compound_match.value.kegg_id,
                             compound_match.key))
        
        return reactants, products
    

class ReactionMatcher(object):
    """Parses reaction queries from users."""
    
    def __init__(self, compound_matcher):
        """Initialize the reaction parser.
        
        Args:
            compound_matcher: a matcher.Matcher object that matches
                              individual compounds.
        """
        self._matcher = compound_matcher
    
    def _MakeReactionCompoundMatch(self, coeff, name):
        compound_matches = self._matcher.Match(name)
        return ReactionCompoundMatch(name, coeff, compound_matches)
    
    def _MakeReactionSide(self, parsed_side):
        side = []
        for coeff, name in parsed_side:
            side.append(self._MakeReactionCompoundMatch(coeff, name))
        return side
    
    def MatchReaction(self, parsed_query):
        """Parse the query as a reaction.
        
        Args:
            parsed_query: query_parser.ParsedReactionQuery object.
        
        Returns:
            An initialized ReactionMatches object.
        """        
        substrates = self._MakeReactionSide(parsed_query.substrates)
        products = self._MakeReactionSide(parsed_query.products)
        
        if not substrates:
            logging.error('Failed to parse substrates.')
            return None
        
        if not products:
            logging.error('Failed to parse products.')
            return None
        
        return ReactionMatches(substrates, products)
        
        

        