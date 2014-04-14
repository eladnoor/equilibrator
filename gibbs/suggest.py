import json
import logging
from django.http import HttpResponse
from django.http import Http404
from gibbs import forms
from gibbs import service_config


def SuggestJson(request):
    """Renders the suggest JSON."""
    form = forms.SearchForm(request.GET)
    if not form.is_valid():
        logging.error(form.errors)
        raise Http404
    
    matcher = service_config.Get().compound_matcher
    query = str(form.cleaned_query)
    matches = matcher.Match(query)
    results = [unicode(m.key) for m in matches]
    types = [m.TypeStr() for m in matches]
    json_data = json.dumps({'query': query,
                            'suggestions': results,
                            'data': types})
    
    return HttpResponse(json_data, mimetype='application/json')