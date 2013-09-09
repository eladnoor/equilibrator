import logging

from gibbs import models
from django.shortcuts import render_to_response


def RefsPage(request):
    """Renders the graph page."""
    sources = models.ValueSource.objects.all()
    sorted_sources = sorted(sources, key=lambda s: s.citation)
    template_data = {"sources": sorted_sources}
    return render_to_response('data_refs.html', template_data)