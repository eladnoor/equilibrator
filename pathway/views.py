import io
import logging
from django.http import HttpResponse
from django.shortcuts import render
from django.template.context_processors import csrf
from .forms import AnalyzePathwayModelForm, BuildPathwayModelForm
from util import constants
from . import pathway_result_page
from . import MaxMinDrivingForce, EnzymeCostMinimization

# PuLP tends to open a TextIOWrapper file without taking care of closing it.
# Therefore, we have to ignore these warnings, but hopefully in future versions
# of PuLP that will be solved.
import warnings
warnings.simplefilter("ignore", ResourceWarning)

def DefinePathwayPage(request):
    """Renders the landing page."""
    template_data = csrf(request)
    
    template_data['ionic_strength'] = \
        request.COOKIES.get('ionic_strength', constants.DEFAULT_IONIC_STRENGTH)
    template_data['ph'] = \
        request.COOKIES.get('pH', constants.DEFAULT_PH)

    response = render(request, 'define_pathway_page.html', template_data)
    return response


def BuildPathwayModel(request):
    """Renders a page for a particular reaction."""
    try:
        form = BuildPathwayModelForm(request.POST, request.FILES)

        if not form.is_valid():
            raise Exception(str(form.errors))

        fp = request.FILES['pathway_file']
        pp, output_fname = pathway_result_page.from_csv(form, fp)
        logging.info(output_fname)

        response = HttpResponse(content_type='text/tab-separated-values')
        response['Content-Disposition'] = 'attachment; filename="%s"' % \
            output_fname
        response.write(pp.to_sbtab())
        
        logging.info('build pathway pH = %.2f' % pp.aq_params.pH)
        response.set_cookie('pH', str(pp.aq_params.pH))
        response.set_cookie('ionic_strength', str(pp.aq_params.ionic_strength))
    
        return response
    except Exception as e:
        logging.error(e)
        template_data = {'pathway': None,
                         'mdf_result': None,
                         'error_message': str(e)}
        return render(request, 'pathway_result_page.html', template_data)


def PathwayResultPage(request):
    """Renders a page for a particular reaction."""
    
    pp = None
    try:
        form = AnalyzePathwayModelForm(request.POST, request.FILES)
        if not form.is_valid():
            raise Exception('Invalid pathway form.')
        
        optimization_method = form.GetOptimizationMethod()
    
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
    
    except ValueError as e:
        logging.error(e)
        template_data = {'pathway': pp,
                         'mdf_result': None,
                         'error_message': str(e)}
        return render(request, 'pathway_result_page.html', template_data)
   
