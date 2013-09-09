#!/usr/bin/env python
# encoding: utf-8

"""
approximate_matching_cmdline.py

Created by Avi Flamholz on 2010-09-21.
Copyright (c) 2010 Weizmann Institute. All rights reserved.
"""

import compound
import approximate_matcher
import reaction_parser
import getopt
import os
import sys


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hoi:v", ["help", "output="])
        except getopt.error, msg:
            raise Usage(msg)
    
        # option processing
        input_filename = "compounds.csv"
        for option, value in opts:
            if option in ("-i", "--infile"):
                input_filename = value
        
        # Read in all the compounds.
        compounds = compound.ReadCompoundsFromCsvFile(input_filename)
        m = approximate_matcher.BackfillingRegexApproxMatcher(compounds,
                                                              max_results=10,
                                                              min_score=0.0)
        rm = reaction_parser.ReactionParser(m)
        
        # REPL loop.
        # TODO(flamholz): Current requires ctrl-C to quit. Allow clean shutdown.
        while True:
            print 'query >',
            
            # Execute the search
            query = raw_input().strip()
            if rm.ShouldParseAsReaction(query):
                print 'Attempting to parse as a reaction.'
                reaction_results = rm.ParseReactionQuery(query)
                if not reaction_results:
                    print 'Failed to parse reaction'
                else:     
                    for substrate in reaction_results.substrates:
                        print substrate.parsed_coeff, substrate.parsed_name,
                        match_names = [m.value.name for m in substrate.matches]
                        print ', '.join(match_names)
                    
                    print '<=>'
                    for product in reaction_results.products:
                        print product.parsed_coeff, product.parsed_name,
                        match_names = [m.value.name for m in product.matches]
                        print ', '.join(match_names)
                    print
            else:  
                print 'Attempting to parse as a single compound.'
                results = m.Match(query.strip())
                
                if not results:
                    print 'Found no results for search: %s' % query
                else:
                    for match in results:
                        print '%s:\t%s\t%f' % (match.value.name,
                                               match.value.kegg_id,
                                               match.score)
                print
                

    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for help use --help"
        return 2


if __name__ == "__main__":
    sys.exit(main())