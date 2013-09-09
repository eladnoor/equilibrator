#!/usr/bin/python

from util import django_utils
import csv
import logging

# NOTE(flamholz): This is crappy. We're using the real database for
# a unit test. I wish I knew of a better way.
django_utils.SetupDjango()

from gibbs import service_config

def main():    
    config = service_config.Get()
    matcher = config.compound_matcher
    
    while True:
        print 'query:',
        query = raw_input()
        matches = matcher.Match(query)
        for m in matches:
            print m.key, m.value.kegg_id, m.score
            

if __name__ == '__main__':
    main()
                
                
