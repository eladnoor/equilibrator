import logging, re

from gibbs import models
from gibbs import constants


class ReactionCompoundMatch(object):
    """A match of a compound in a reaction.
    
    Contains the parsed name (what the user entered), the parsed
    stoichiometric coefficient, and a list of matcher.Match objects
    of the potential matches.
    """
    def __init__(self, parsed_name, parsed_coeff, parsed_phase, matches):
        self.parsed_name = parsed_name
        self.parsed_coeff = parsed_coeff
        self.parsed_phase = parsed_phase
        self.matches = matches
    
    def __str__(self):
        return '%d %s(%s), matches: %s' % (self.parsed_coeff,
                                           self.parsed_name,
                                           self.parsed_phase,
                                           ', '.join(self.matches))
    
    def ParsedDataEqual(self, other):
        """Checks if the parsed data (name, coefficient) are equal
           for two ReactionCompoundMathes.
        """
        return (self.parsed_coeff == other.parsed_coeff and
                self.parsed_name == other.parsed_name and
                self.parsed_phase == other.parsed_phase)

    
class ReactionMatches(object):
    """A reaction parsed from a query with all possible matches."""
    
    def __init__(self, reactants=None):
        """Initialize the ReactionMatches object.
        
        Args:
            reactants: a list of ReactionCompoundMatches for the reactants.
        """
        self.reactants = reactants or []
    
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
        for c in self.reactants:
            compound_match = self._FindFirstCompoundMatch(c.matches)
            if not compound_match:
                return None
            
            reactants.append({'coeff': c.parsed_coeff,
                              'kegg_id': compound_match.value.kegg_id,
                              'phase': c.parsed_phase,
                              'name': compound_match.key})
        
        return reactants
    

class ReactionMatcher(object):
    """Parses reaction queries from users."""
    
    def __init__(self, compound_matcher):
        """Initialize the ReactionMatcher.
        
        Args:
            compound_matcher: a matcher.Matcher object that matches
                individual compounds.
        """
        self._matcher = compound_matcher
    
    def _MakeReactionCompoundMatch(self, coeff, name):
        compound_name, phase_suffix = ReactionMatcher._StripPhases(name)
        logging.debug("Name = %s, phase = %s" % (compound_name, phase_suffix))
        compound_matches = self._matcher.Match(compound_name)
        return ReactionCompoundMatch(compound_name, coeff, phase_suffix,
                                     compound_matches)

    @staticmethod
    def _StripPhases(name):
        compound_name = name
        phase_name = None

        m = re.match('(.*)\((aq|s|l|g)\)$', name)
        if m is not None:
            compound_name, phase_subscript = m.groups()
            phase_name = constants.PHASE_SUBSCRIPT_TO_NAME[phase_subscript]

        return compound_name, phase_name
    
    def MatchReaction(self, parsed_query):
        """Parse the query as a reaction.
        
        Args:
            parsed_query: query_parser.ParsedReactionQuery object.
        
        Returns:
            An initialized ReactionMatches object.
        """  
        reactants = []
        for coeff, name in parsed_query.substrates:
            reactants.append(self._MakeReactionCompoundMatch(-coeff, name))
        for coeff, name in parsed_query.products:
            reactants.append(self._MakeReactionCompoundMatch(coeff, name))
        
        if not reactants:
            logging.error('Failed to parse reaction.')
            return None
        
        return ReactionMatches(reactants)
        
        

        