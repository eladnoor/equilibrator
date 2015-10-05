from django.shortcuts import render_to_response
from django.core.context_processors import csrf
from gibbs import constants

def DefinePathwayPage(request):
    """Renders the landing page."""
    template_params = csrf(request)
    template_params['constants'] = constants
    return render_to_response('define_pathway_page.html', template_params)