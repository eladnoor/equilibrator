import logging
import json

from django.http import HttpResponse
from django.http import HttpResponseBadRequest
from django.http import HttpResponseServerError
from django.views.decorators.csrf import csrf_exempt
from gibbs import constants
from gibbs import models
from gibbs import reaction


MAX_REACTIONS = 50

@csrf_exempt
def ReactionAPI(request):    
    """Outputs JSON data for up to MAX_REACTIONS reactions."""
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
        return HttpResponseBadRequest('Request is not valid JSON.')
    
    try:
        reaction_ids = parsed_data.get('KEGG_reactions', [])
        reaction_ids = map(str, reaction_ids)
    except Exception, e:
        logging.warning(e)
        HttpResponseBadRequest('Request includes invalid KEGG reaction IDs.')
        
    if len(reaction_ids) > MAX_REACTIONS:
        logging.info('API request too large, ignoring.')
        return HttpResponseBadRequest('Requested more than %d reactions.' %
                                      MAX_REACTIONS)
    
    ph = parsed_data.get('pH', constants.DEFAULT_PH)
    i_s = parsed_data.get('ionic_strength', constants.DEFAULT_IONIC_STRENGTH)
    
    # Fetch reactions from DB.
    try:
        stored_reactions = []
        if reaction_ids:
            stored_reactions = models.StoredReaction.objects.select_related(
                ).filter(kegg_id__in=reaction_ids)
    except Exception, e:
        # TODO(flamholz): catch more specific errors
        logging.error(e)
        return HttpResponseServerError('DB Query failed. Try again.')
        
    
    # Make them into JSON!
    reactions = [reaction.Reaction.FromStoredReaction(sr, pH=ph, ionic_strength=i_s)
                 for sr in stored_reactions]
    
    json_rxns = [r.ToJson() for r in reactions]
    json_dict = {'reactions': json_rxns,
                 'pH': ph,
                 'ionic_strength': i_s}
    json_data = json.dumps(json_dict)
    return HttpResponse(json_data, mimetype='application/json')