import logging, io

from django.http import HttpResponseBadRequest
from django.shortcuts import render_to_response
from django.views.decorators.csrf import csrf_exempt
from gibbs.forms import PathwayForm
from pathways import ParsedPathway


def PathwayResultPage(request):    
    """Renders a page for a particular reaction."""
    form = PathwayForm(request.POST, request.FILES)
    if not form.is_valid():
        logging.error(form.errors)
        return HttpResponseBadRequest('Invalid pathway form.')
    
    f_data = unicode(request.FILES['pathway_file'].read())
    sio = io.StringIO(f_data, newline=None)  # universal newline mode
    path = ParsedPathway.from_file(sio)
    mdf_result = path.calc_mdf()
    template_data = {'pathway': path,
                     'mdf_result': mdf_result}

    return render_to_response('pathway_result_page.html', template_data)
    