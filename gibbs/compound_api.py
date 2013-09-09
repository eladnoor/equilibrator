import logging
import json

from django.http import HttpResponse
from django.http import HttpResponseBadRequest
from django.http import HttpResponseServerError
from django.views.decorators.csrf import csrf_exempt
from gibbs import constants
from gibbs import models
from gibbs import reaction


MAX_COMPOUNDS = 50

@csrf_exempt
def CompoundAPI(request):
    """Outputs JSON data for up to MAX_COMPOUNDS compounds."""
    data = request.raw_post_data
    if not data:
        return HttpResponseBadRequest('No request data.')
    
    try:
        parsed_data = json.loads(data)
    except Exception, e:
        logging.warning(e)
        return HttpResponseBadRequest('Request is not valid JSON.')
    
    if not parsed_data:
        logging.warning('Empty API query.')
        return HttpResponseBadRequest('Request is not valid.')
    
    try:
        inchis = map(str, parsed_data.get('InChI_identifiers', []))
        kegg_ids = map(str, parsed_data.get('KEGG_compounds', []))
    except Exception, e:
        logging.warning(e)
        HttpResponseBadRequest('Request includes invalid InChI compound IDs.')
    
    if len(inchis) + len(kegg_ids) > MAX_COMPOUNDS:
        logging.info('API request too large, ignoring.')
        return HttpResponseBadRequest('Requested more than %d reactions.' %
                                      MAX_COMPOUNDS)
    
    ph = parsed_data.get('pH', constants.DEFAULT_PH)
    i_s = parsed_data.get('ionic_strength', constants.DEFAULT_IONIC_STRENGTH)
    
    # Fetch compounds from DB.
    stored_compounds = []
    try:
        if inchis:
            stored_compounds.extend(models.Compound.objects.select_related()
                                    .filter(inchi__in=inchis))
        if kegg_ids:
            stored_compounds.extend(models.Compound.objects.select_related()
                                    .filter(kegg_id__in=kegg_ids))
    except Exception, e:
        # TODO(flamholz): catch more specific errors
        logging.error(e)
        return HttpResponseServerError('DB Query failed. Try again.')
        
    json_compounds = [c.ToJson() for c in stored_compounds]
    json_dict = {'compounds': json_compounds,
                 'pH': ph,
                 'ionic_strength': i_s}
    json_data = json.dumps(json_dict)
    return HttpResponse(json_data, mimetype='application/json')
