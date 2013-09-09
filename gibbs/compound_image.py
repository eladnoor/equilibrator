import logging
import json

from django.http import HttpResponse
from django.http import HttpResponseBadRequest
from django.http import HttpResponseServerError
from django.views.decorators.csrf import csrf_exempt
from gibbs import constants
from gibbs import models
from gibbs import reaction
import base64
import os
from equilibrator.settings import MEDIA_ROOT
from django.core.files import File

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
    if not compound.thumbnail:
        image_data = open(NO_STRUCTURE_THUMBNAIL, 'r').read()
    else:
        image_data = base64.decodestring(compound.thumbnail)
    return HttpResponse(image_data, mimetype='image/png')
