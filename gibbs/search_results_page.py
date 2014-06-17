import logging

from django.http import Http404
from django.shortcuts import render_to_response
from gibbs import conditions
from gibbs import reaction
from gibbs.forms import SearchForm
from gibbs import service_config
from profiling import profile

@profile    
def ResultsPage(request):
    """Renders the search results page for a given query."""
    form = SearchForm(request.GET)
    if not form.is_valid():
        raise Http404
    
    logging.debug('Generating a search result page')
    query_parser = service_config.Get().query_parser
    reaction_matcher = service_config.Get().reaction_matcher
    matcher = service_config.Get().compound_matcher
    
    query = form.cleaned_query
    if not query.strip():
        response = render_to_response('main.html', {})
        return response
    
    # Check if we should parse and process the input as a reaction.
    if query_parser.IsReactionQuery(query):
        logging.debug('Parsing the query as a reaction')
        try:
            parsed_reaction = query_parser.ParseReactionQuery(query)
        except Exception:
            return render_to_response('parse_error_page.html')

        reaction_matches = reaction_matcher.MatchReaction(parsed_reaction)
        best_reaction = reaction_matches.GetBestMatch()
        
        if not best_reaction:
            return render_to_response('search_error_page.html')

        logging.debug('Generating a reaction from the matched KEGG IDs')
        rxn = reaction.Reaction.FromIds(best_reaction)
        rxn.aq_params = conditions.AqueousParams.FromForm(form, request.COOKIES) 
        logging.info('Aqueous parameters: ' + str(rxn.aq_params))
        
        response = render_to_response('reaction_page.html', rxn.GetTemplateData(query))
        return response

    else:
        # Otherwise we try to parse it as a single compound.
        logging.debug('Parsing the query as a single compound')
        results = matcher.Match(query)
        template_data = {}        
        template_data['compound_results'] = [m for m in results if m.IsCompound()]
        template_data['enzyme_results'] = [m for m in results if m.IsEnzyme()]
        template_data['enzymes_first'] = results and results[0].IsEnzyme()
            
        response = render_to_response('search_results.html', template_data)
        return response

    raise Http404
