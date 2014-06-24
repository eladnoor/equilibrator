from django.shortcuts import render_to_response
from equilibrator.gibbs import constants
from django.template.context import RequestContext


def AboutPage(request):
    """Renders the about page."""
    return render_to_response('about.html', {})


def FAQPage(request):
    """Renders the FAQ page."""
    return render_to_response('faq.html', {})

def WhatsNewPage(request):
    """Renders the FAQ page."""
    return render_to_response('new_in_2_0.html', {})


def CitePage(request):
    """Renders the FAQ page."""
    return render_to_response('cite.html', {})


def RedoxReview(request):
    """Renders robots.txt."""
    return render_to_response('redox_review.html', {})


def ClassicReactions(request):
    """Renders the classic reactions page."""
    return render_to_response('classic_reactions.html', {})


def DownloadPage(request):
    """Renders the download page."""
    
    ph_values = map(lambda x: '%.1f' % x, constants.PH_RANGE_VALUES)
    return render_to_response('download.html', {'ph_values': ph_values},
                              context_instance = RequestContext(request))
    

def Robots(request):
    """Renders robots.txt."""
    return render_to_response('robots.txt', {})