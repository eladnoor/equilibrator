from matching import approximate_matcher
from matching import query_parser
from matching import reaction_matcher
from util import singleton


@singleton.Singleton
class ServiceConfig(object):
    """
        A singleton class that contains global service configuration and state.

        All state/config is considered immutable.
            query_parser: parses search queries.
            search_matcher: provides search results and completions.
            reaction_matcher: matches reactions against the database.
    """
    def __init__(self):
        self._query_parser = query_parser.QueryParser()
        self._search_matcher = approximate_matcher.CascadingMatcher(
            max_results=10, min_score=0.1,
            match_enzymes=True, return_fast=False)

        # Don't want to match enzymes when matching reactions.
        self._single_compound_matcher = approximate_matcher.CascadingMatcher(
            max_results=1, min_score=0.1,
            match_enzymes=False, return_fast=True)
        self._reaction_matcher = reaction_matcher.ReactionMatcher(
            self._single_compound_matcher)

    query_parser = property(lambda self: self._query_parser)
    search_matcher = property(lambda self: self._search_matcher)
    reaction_matcher = property(lambda self: self._reaction_matcher)


def Get():
    """
        Convenience method to get the single ServiceConfig instance.
    """
    return ServiceConfig()  # Singleton decorator ensures there's only one
