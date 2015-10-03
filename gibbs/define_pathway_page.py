from django.shortcuts import render_to_response
from django.core.context_processors import csrf

def DefinePathwayPage(request):
    """Renders the landing page."""
    template_params = csrf(request)
    return render_to_response('define_pathway_page.html', template_params)