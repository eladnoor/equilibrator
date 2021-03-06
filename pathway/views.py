import io
import logging
import os
from django.http import HttpResponse, HttpResponseBadRequest
from django.shortcuts import render
from django.template.context_processors import csrf
from .forms import AnalyzePathwayModelForm, BuildPathwayModelForm
from util import constants
from . import pathway_result_page
from . import ParsedPathway, PathwayParseError

def DefinePathwayPage(request):
    """Renders the landing page."""
    template_params = csrf(request)
    template_params['constants'] = constants
    return render(request, 'define_pathway_page.html', template_params)


def PathwayResultPage(request):
    """Renders a page for a particular reaction."""
    form = AnalyzePathwayModelForm(request.POST, request.FILES)
    if not form.is_valid():
        logging.error(form.errors)
        return HttpResponseBadRequest('Invalid pathway form.')

    try:
        f_data = str(request.FILES['pathway_file'].read(), encoding="ascii")
        sio = io.StringIO(f_data, newline=None)  # universal newline mode
        reactions, fluxes, keqs, bounds = pathway_result_page.read_sbtabs(sio)
        pp = ParsedPathway.from_full_sbtab(
            reactions, fluxes, bounds, keqs)
        logging.info('Parsed pathway.')
    except PathwayParseError as ppe:
        logging.error(ppe)
        return HttpResponseBadRequest(ppe.message)
    except Exception as e:
        logging.error(e)
        template_data = {'pathway': None,
                         'mdf_result': None,
                         'error_message': str(e)}
        return render(request, 'pathway_result_page.html', template_data)

    if len(pp.reactions) == 0:
        logging.error('Pathway contains no reactions')
        template_data = {'pathway': pp,
                         'mdf_result': None,
                         'error_message': 'Empty pathway'}
        return render(request, 'pathway_result_page.html', template_data)

    try:
        # calculate the MDF with the specified bounds. Render template.
        mdf_result = pp.calc_mdf()
        template_data = {'pathway': pp,
                         'mdf_result': mdf_result}
        logging.info('Calculated MDF %s', mdf_result.mdf)
        return render(request, 'pathway_result_page.html', template_data)
    except Exception as e:
        logging.error(e)
        template_data = {'pathway': pp,
                         'mdf_result': None,
                         'error_message': str(e)}
        return render(request, 'pathway_result_page.html', template_data)

def BuildPathwayModel(request):
    """Renders a page for a particular reaction."""
    form = BuildPathwayModelForm(request.POST, request.FILES)
    if not form.is_valid():
        logging.error(form.errors)
        return HttpResponseBadRequest('Invalid pathway form.')

    try:
        bounds = pathway_result_page.make_bounds(form)
        aq_params = pathway_result_page.make_aq_params(form)
    except Exception as e:
        logging.error(e)
        return HttpResponseBadRequest(e)

    try:
        f = request.FILES['pathway_file']
        fname_base, ext = os.path.splitext(f.name)
        output_fname = '%s_pH%.2f_I%.2f.tsv' % (
            fname_base, aq_params.pH, aq_params.ionic_strength)
        logging.info(output_fname)

        pp = ParsedPathway.from_csv_file(
            f, bounds=bounds, aq_params=aq_params)
    except PathwayParseError as ppe:
        logging.error(ppe)
        return HttpResponseBadRequest(ppe)

    response = HttpResponse(content_type='text/tab-separated-values')
    response['Content-Disposition'] = 'attachment; filename="%s"' % \
        output_fname
    response.write(pp.to_full_sbtab())

    return response
