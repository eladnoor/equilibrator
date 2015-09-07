import json
import logging
from django.http import HttpResponse
from django.http import Http404
from gibbs import forms
from gibbs import service_config
from matching.approximate_matcher import HaystackApproxMatcher
from haystack.query import SearchQuerySet


def SuggestJson(request):
    """Renders the suggest JSON."""
    form = forms.SuggestForm(request.GET)
    if not form.is_valid():
        logging.error(form.errors)
        raise Http404
      
    query = str(form.cleaned_query)
    suggestions = []
    if query:
        matcher = service_config.Get().compound_matcher
        matches = matcher.Match(query)
        suggestions = [{'value': str(m.key), 'data': {'cat': m.key.TypeStr()}}
                       for m in matches]
    output = {'query': query,
              'suggestions': suggestions}
    json_data = json.dumps(output)
    
    return HttpResponse(json_data, mimetype='application/json')