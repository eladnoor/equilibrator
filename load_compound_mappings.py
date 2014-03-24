import json
import logging
import gzip

import util.django_utils

from gibbs import models

DEFAULT_COMPOUND_MAPPING_FILENAME = 'data/compound_replacement_mapping.json.gz'

def LoadEquivalentCompounds(json_filename=DEFAULT_COMPOUND_MAPPING_FILENAME):
    """Loads mappings between compounds."""
    parsed_json = json.load(gzip.open(json_filename))
    
    for mapping in parsed_json:
        from_kegg_id = mapping['FROM ID']
        to_kegg_id = mapping['TO ID']
        
        try:
            from_compound = models.Compound.objects.get(kegg_id=from_kegg_id)
            to_compound = models.Compound.objects.get(kegg_id=to_kegg_id)
            from_compound.replace_with = to_compound
            from_compound.save()
        except Exception, e:
            logging.error(e)
            continue


if __name__ == '__main__':
    LoadEquivalentCompounds()