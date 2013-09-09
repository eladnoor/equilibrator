# -*- coding: utf-8 -*-

import logging
import pyparsing
import re


class ParseError(Exception):
    pass


def _parsedCompound(c_list):
    """Always put a stoichiometric coefficient with a compound."""
    if len(c_list) == 2:
        return c_list[0], c_list[1]
    return 1, c_list[0]


def _MakeReactionParser():
    """Builds a pyparsing-based recursive descent parser for chemical reactions."""
    
    # Coefficients are usually integral, but they can be floats or fractions too.
    int_coeff = pyparsing.Word(pyparsing.nums)
    float_coeff = pyparsing.Word(pyparsing.nums + '.' + pyparsing.nums)
    frac_coeff = int_coeff + '/' + int_coeff
    int_coeff.setParseAction(lambda i:int(i[0]))
    float_coeff.setParseAction(lambda t:float(t[0]))
    frac_coeff.setParseAction(lambda f:float(f[0])/float(f[2]))
    
    coeff = pyparsing.Or([int_coeff, float_coeff, frac_coeff])
    optional_coeff = pyparsing.Optional(coeff)
    
    compound_separator = pyparsing.Literal('+').suppress()
    
    compound_name_component = pyparsing.Word(pyparsing.alphanums + "()",
                                             pyparsing.alphanums + "-+,()'")
    compound_name = pyparsing.Forward()
    compound_name << (compound_name_component + pyparsing.ZeroOrMore(compound_name_component))
    compound_name.setParseAction(lambda s: ' '.join(s))
    
    compound_with_coeff = pyparsing.Forward()
    compound_with_coeff << ((optional_coeff + compound_name) | compound_name)
    compound_with_coeff.setParseAction(_parsedCompound)
    compound_with_coeff.setResultsName("compound")
    
    compound_with_separator = pyparsing.Forward()
    compound_with_separator << (compound_with_coeff + compound_separator)
    
    reaction_side = pyparsing.Forward()
    reaction_side << (pyparsing.ZeroOrMore(compound_with_separator) +
                      compound_with_coeff)
    reaction_side.setParseAction(lambda l: [l])
    reaction_side.setResultsName("reaction_side")
    
    side_separators = [pyparsing.Literal(s) for s in ("=", "->", "=>", "<=>", u'→')]
    side_separator = pyparsing.Or(side_separators).suppress()
    
    reaction = pyparsing.Forward()
    reaction << (reaction_side + side_separator + reaction_side)
    return reaction


class ParsedReactionQuery(object):
    """A parsed reaction query."""
    
    def __init__(self, substrates=None, products=None):
        """Initialize the ParsedReaction object.
        
        Args:
            reactants: a list of tuples for the reactants.
            products: a list of tuples for the products.
        """
        self.substrates = substrates or []
        self.products = products or []
    
    def __eq__(self, other):
        """Equality test."""
        r = frozenset(self.substrates)
        p = frozenset(self.products)
        o_r = frozenset(other.substrates)
        o_p = frozenset(other.products)
        
        reactants_diff = r.symmetric_difference(o_r)
        products_diff = p.symmetric_difference(o_p)
        
        if not reactants_diff and not products_diff:
            return True
        
        return False
    
    def __str__(self):
        joined_rs = ['%s %s' % (c,r) for c,r in self.substrates]
        joined_ps = ['%s %s' % (c,p) for c,p in self.products]
        return '%s => %s' % (' + '.join(joined_rs), ' + '.join(joined_ps))
    

class QueryParser(object):
    """Parses search queries."""
    
    REACTION_PATTERN = u'.*(=>|<=>|=|->|<->|→).*'
    REACTION_MATCHER = re.compile(REACTION_PATTERN)
    
    def __init__(self):
        """Initialize the parser."""
        self._rparser = _MakeReactionParser()
        
    def IsReactionQuery(self, query):
        """Returns True if this query is likely to be a reaction query.
        
        Args:
            query: the query string.
        """
        m = self.REACTION_MATCHER.match(query.strip())
        return m is not None
    
    def ParseReactionQuery(self, query):
        """Parse the query as a reaction.
        
        Args:
            query: the query string.
        
        Returns:
            An initialized ParsedReaction object, or None if parsing failed.
        """
        try:
            results = self._rparser.parseString(query)
            substrates, products = results
            return ParsedReactionQuery(substrates, products)
        except pyparsing.ParseException,msg:
            logging.error('Failed to parse query %s', query)
            raise ParseError(msg)
        