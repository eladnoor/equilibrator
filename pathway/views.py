import io
import logging
import os
from django.http import HttpResponse, HttpResponseBadRequest
from django.shortcuts import render
from django.template.context_processors import csrf
from .forms import AnalyzePathwayModelForm, BuildPathwayModelForm
from util import constants
from . import pathway_result_page
from . import MaxMinDrivingForce, EnzymeCostMinimization, PathwayParseError

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

    pp = None
    try:
        optimization_method = form.cleaned_data['optimization_method']
        assert optimization_method in ['MDF', 'ECM']
    
        f_data = str(request.FILES['pathway_file'].read(), encoding="ascii")
        sio = io.StringIO(f_data, newline=None)  # universal newline mode
        if optimization_method == 'MDF':
            pp = MaxMinDrivingForce.from_sbtab_file(sio)
        else:
            pp = EnzymeCostMinimization.from_sbtab_file(sio)

        logging.info('Parsed pathway.')
        
        if not pp.validate_dGs():
            raise Exception('Supplied reaction dG values are inconsistent '
                            'with the stoichiometric matrix.')

        if pp.is_empty():
            raise Exception('Empty pathway')

        # calculate the MDF with the specified bounds. Render template.
        path_data = pp.analyze()
        logging.info('Analyzed pathway.')
        template_data = {'pathway': pp,
                         'path_data': path_data}
        logging.info('Calculated %s score: %s' %
                     (optimization_method, path_data.score))
        return render(request, 'pathway_result_page.html', template_data)
    
    except IOError as e:
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
        optimization_method = form.cleaned_data['optimization_method']
        assert optimization_method in ['MDF', 'ECM']
    except Exception as e:
        logging.error(e)
        return HttpResponseBadRequest(e)

    try:
        fp = request.FILES['pathway_file']
        fname_base, ext = os.path.splitext(fp.name)
        output_fname = '%s_pH%.2f_I%.2f_%s.tsv' % (
            fname_base, aq_params.pH, aq_params.ionic_strength, optimization_method)
        logging.info(output_fname)

        if optimization_method == 'MDF':
            pp = MaxMinDrivingForce.from_csv(
                fp, bounds=bounds, aq_params=aq_params)
        else:
            pp = EnzymeCostMinimization.from_csv(
                fp, bounds=bounds, aq_params=aq_params)
    except PathwayParseError as ppe:
        logging.error(ppe)
        return HttpResponseBadRequest(ppe)

    response = HttpResponse(content_type='text/tab-separated-values')
    response['Content-Disposition'] = 'attachment; filename="%s"' % \
        output_fname
    response.write(pp.to_sbtab())

    return response
