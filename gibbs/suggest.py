import json
import logging
from django.http import HttpResponse
from django.http import Http404
from gibbs import forms
from gibbs import service_config
from haystack.query import SearchQuerySet


def SuggestJson(request):
    """Renders the suggest JSON."""
    form = forms.SearchForm(request.GET)
    if not form.is_valid():
        logging.error(form.errors)
        raise Http404
      
    query = str(form.cleaned_query)  
    sqs = SearchQuerySet().autocomplete(name_auto=query)[:10]
    suggestions = [result.text for result in sqs]
    common_names = [result.object for result in sqs]
    types = [cn.TypeStr() for cn in common_names]
    json_data = json.dumps({'query': query,
                            'suggestions': suggestions,
                            'data': types})
    
    return HttpResponse(json_data, mimetype='application/json')