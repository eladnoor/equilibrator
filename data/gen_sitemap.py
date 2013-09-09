#!/usr/bin/python

from util import django_utils

# NOTE(flamholz): This is crappy. We're using the real database for
# a unit test. I wish I knew of a better way.
django_utils.SetupDjango()

from gibbs import models

BASE_URL = 'http://equilibrator.weizmann.ac.il'
LINK_FORMAT = BASE_URL + '%s\n'

def main():
    urls = [BASE_URL + '/\n',
            BASE_URL + '/faq\n',
            BASE_URL + '/about\n',
            BASE_URL + '/classic_reactions\n',
            BASE_URL + '/download\n',
            BASE_URL + '/data_refs\n']
    
    for compound in models.Compound.objects.all():
        urls.append(LINK_FORMAT % compound.link)
        
    for reaction in models.StoredReaction.objects.all():
        urls.append(LINK_FORMAT % reaction.link)
    
    for enzyme in models.Enzyme.objects.all():
        urls.append(LINK_FORMAT % enzyme.link)
    
    
    f = open('sitemap.txt', 'w')
    f.writelines(urls)
    f.close()
            

if __name__ == '__main__':
    main()
                
                
