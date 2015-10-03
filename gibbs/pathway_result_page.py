import logging
import io

from django.http import HttpResponseBadRequest
from django.shortcuts import render_to_response
from django.views.decorators.csrf import csrf_exempt
from gibbs.forms import PathwayForm
from os import path
from pathways import ParsedPathway
from pathways.bounds import Bounds

RELPATH = path.dirname(path.realpath(__file__))
COFACTORS_FNAME = path.join(RELPATH, '../pathways/data/cofactors.csv')


def make_bounds(request, form):
    max_c = form.cleaned_data.get('max_c')
    min_c = form.cleaned_data.get('min_c')
    default_ub = max_c or 100e-3
    default_lb = min_c or 1e-6
    conc_f = request.FILES.get('bounds_file')
    logging.info('default_lb = %.2g, default_ub = %.2g.',
                 default_lb, default_ub)

    if conc_f:
        logging.info('Reading concentrations from uploaded file.')
        concs_data = unicode(conc_f.read())
        sio = io.StringIO(concs_data, newline=None)  # universal newline mode
        bounds = Bounds.from_csv_file(
            sio, default_lb=default_lb, default_ub=default_ub)
        return bounds
    if max_c or min_c:
        logging.info('Reading default concentration bounds.')
        bounds = Bounds.from_csv_filename(
            COFACTORS_FNAME, default_lb=default_lb, default_ub=default_ub)
        return bounds


def PathwayResultPage(request):    
    """Renders a page for a particular reaction."""
    form = PathwayForm(request.POST, request.FILES)
    if not form.is_valid():
        logging.error(form.errors)
        return HttpResponseBadRequest('Invalid pathway form.')

    f_data = unicode(request.FILES['pathway_file'].read())
    sio = io.StringIO(f_data, newline=None)  # universal newline mode
    path = ParsedPathway.from_file(sio)

    # If specific bounds are supplied, use them.
    bounds = make_bounds(request, form)
    path.set_bounds(bounds)

    # calculate the MDF with the specified bounds. Render template.
    mdf_result = path.calc_mdf()
    template_data = {'pathway': path,
                     'mdf_result': mdf_result}
    return render_to_response('pathway_result_page.html', template_data)
    