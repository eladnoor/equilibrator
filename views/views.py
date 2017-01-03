import logging
import json
from django.http import HttpResponse
from django.shortcuts import render
from django.http import Http404
from django.apps import apps
from gibbs import conditions, constants, service_config
from gibbs.forms import EnzymeForm, SearchForm, SuggestForm
import util.django_utils


def DownloadPage(request):
    """Renders the download page."""
    ph_values = map(lambda x: '%.1f' % x, constants.PH_RANGE_VALUES)
    return render(request, 'download.html', {'ph_values': ph_values})


def RefsPage(request):
    """Renders a list of relevant literature."""
    sources = apps.get_model('gibbs.ValueSource').objects.all()
    sorted_sources = sorted(sources, key=lambda s: s.citation)
    template_data = {"sources": sorted_sources}
    return render(request, 'data_refs.html', template_data)


def EnzymePage(request):
    """Renders a page for a particular enzyme."""
    form = EnzymeForm(request.GET)
    if not form.is_valid():
        logging.error(form.errors)
        raise Http404

    enz = apps.get_model('gibbs.Enzyme').objects.get(ec=form.cleaned_ec)
    template_data = {'is_superuser': util.django_utils.IsSuperUser(request),
                     'enzyme': enz}
    return render(request, 'enzyme_page.html', template_data)


def ResultsPage(request):
    """Renders the search results page for a given query."""
    form = SearchForm(request.GET)
    if not form.is_valid():
        raise Http404
    logging.debug('Generating a search result page')
    query_parser = service_config.Get().query_parser
    reaction_matcher = service_config.Get().reaction_matcher
    matcher = service_config.Get().search_matcher
    query = form.cleaned_query
    logging.debug('Query: "%s"', query)
    if not query.strip():
        response = render(request, 'main.html', {})
        return response
    # Check if we should parse and process the input as a reaction.
    if query_parser.IsReactionQuery(query):
        logging.info('Parsing the query as a reaction')
        try:
            parsed_reaction = query_parser.ParseReactionQuery(query)
        except Exception:
            return render(request, 'parse_error_page.html')

        reaction_matches = reaction_matcher.MatchReaction(parsed_reaction)
        best_reaction = reaction_matches.GetBestMatch()
        if not best_reaction:
            return render(request, 'search_error_page.html', {'query': query})

        logging.debug('Generating a reaction from the matched KEGG IDs')
        aq_params = conditions.AqueousParams.FromForm(form, request.COOKIES)
        rxn = apps.get_model('gibbs.Reaction').FromIds(best_reaction, aq_params)
        response = render(request, 'reaction_page.html',
                          rxn.GetTemplateData(query))
        return response

    else:
        # Otherwise we try to parse it as a single compound.
        logging.debug('Parsing the query as a single compound/enzyme')
        results = matcher.Match(query)
        template_data = {}
        template_data['compound_results'] = \
            [m for m in results if m.IsCompound()]
        template_data['enzyme_results'] = [m for m in results if m.IsEnzyme()]
        template_data['enzymes_first'] = results and results[0].IsEnzyme()
        template_data['query'] = query
        response = render(request, 'search_results.html', template_data)
        return response

    raise Http404


def SuggestJson(request):
    """Renders the suggest JSON."""
    form = SuggestForm(request.GET)
    if not form.is_valid():
        logging.error(form.errors)
        raise Http404
    query = str(form.cleaned_query)
    suggestions = []
    if query:
        matcher = service_config.Get().search_matcher
        matches = matcher.Match(query)
        suggestions = [{'value': str(m.key), 'data': {'cat': m.key.TypeStr()}}
                       for m in matches]
    output = {'query': query,
              'suggestions': suggestions}
    json_data = json.dumps(output)
    return HttpResponse(json_data, content_type='application/json')
