import logging
import json

from django.http import HttpResponse
from django.http import HttpResponseBadRequest
from django.views.decorators.csrf import csrf_exempt
import types
from equilibrator.util.inchi_api_test import CreateJSON
from toolbox.database import SqliteDatabase
from pygibbs.unified_group_contribution_util import InChI2FormationEnergy
from pygibbs.group_decomposition import GroupDecompositionError

MAX_COMPOUNDS = 50

@csrf_exempt
def InChIAPI(request):
    """Outputs JSON data for up to MAX_COMPOUNDS compounds."""
    data = request.raw_post_data
    if not data:
        #return HttpResponseBadRequest('No request data.')
        data = CreateJSON()
    
    try:
        parsed_data = json.loads(data)
    except Exception, e:
        logging.warning(e)
        return HttpResponseBadRequest('Request is not valid JSON.')
    
    if not parsed_data or not type(parsed_data) == types.ListType:
        logging.warning('Empty API query.')
        return HttpResponseBadRequest('Request is not valid.')
    
    queue = []
    for d in parsed_data:
        if 'name' in d and 'InChI' in d and 'pKas' in d:
            name = str(d['name'])
            inchi = str(d['InChI'])
            try:
                dissociation_constants = [float(x) for x in d['pKas']]
                queue.append((name, inchi, dissociation_constants))
            except Exception, e:
                logging.warning(e)
                HttpResponseBadRequest('Request includes invalid pKa values')
    
    if len(queue) > MAX_COMPOUNDS:
        logging.info('API request too large, ignoring.')
        return HttpResponseBadRequest('Requested more than %d reactions.' %
                                      MAX_COMPOUNDS)
    
    result = []

    inchi2dg = InChI2FormationEnergy()
    db_out = SqliteDatabase('/home/eladn/workspace/milo-lab/src/equilibrator/data/ugc.sqlite', 'w')
    inchi2dg.FromDatabase(db_out)

    for name, inchi, pKas in queue:
        d = {'name': name, 'InChI': inchi, 'pKas': pKas}
        try:
            dG0, nH, charge, nMg, ker = inchi2dg.EstimateInChI(inchi)
            ker = ker.round(10)
            d['kernel'] = inchi2dg.ArrayToSparseRep(ker)
            d['pseudoisomers'] = inchi2dg.GenerateAllPseudoisomers(dG0, nH,
                                                            charge, nMg, pKas)
        except GroupDecompositionError:
            d['error'] = "Cannot decompose this compound into groups"
        
        result.append(d)
        
    json_data = json.dumps(result, indent=4)
    return HttpResponse(json_data, mimetype='application/json')
