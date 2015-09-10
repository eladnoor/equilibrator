from django.http import HttpResponse
from django.http import HttpResponseBadRequest
from django.views.decorators.csrf import csrf_exempt
from gibbs import models
import base64
from equilibrator.settings import MEDIA_ROOT

NO_STRUCTURE_THUMBNAIL = '/'.join([MEDIA_ROOT, 'images', 'structure_not_available.png'])

@csrf_exempt
def CompoundImage(request):
    data = request.GET
    compound_id = data.get('compoundId', None)
    if compound_id is None:
        return HttpResponseBadRequest('No request data.')
    
    compounds = models.Compound.objects.filter(kegg_id=compound_id)
    if not compounds:
        return HttpResponseBadRequest('No such compound.')
    
    compound = compounds[0]
    if not compound.thumbnail or compound.thumbnail == 'error':
        image_data = open(NO_STRUCTURE_THUMBNAIL, 'r').read()
    else:
        image_data = base64.decodestring(compound.thumbnail)
    return HttpResponse(image_data, content_type='image/png')
