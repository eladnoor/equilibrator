import logging

from django.http import Http404
from django.shortcuts import render_to_response
from gibbs import conditions
from gibbs import reaction
from gibbs import search_form
from gibbs import service_config

    
def ResultsPage(request):
    """Renders the search results page for a given query."""
    form = search_form.SearchForm(request.GET)
    if not form.is_valid():
        raise Http404
    
    logging.debug('Generating a search result page')
    query_parser = service_config.Get().query_parser
    reaction_matcher = service_config.Get().reaction_matcher
    matcher = service_config.Get().compound_matcher
    
    query = form.cleaned_query
    ph = form.cleaned_ph
    pmg = form.cleaned_pmg
    ionic_strength = form.cleaned_ionic_strength
    template_data = {'query': query,
                     'ph': ph,
                     'pmg': pmg, 
                     'ionic_strength': ionic_strength}
        
    # Check if we should parse and process the input as a reaction.
    if query_parser.IsReactionQuery(query):
        logging.debug('Parsing the query as a reaction')
        try:
            parsed_reaction = query_parser.ParseReactionQuery(query)
        except Exception:
            return render_to_response('parse_error_page.html', template_data)

        reaction_matches = reaction_matcher.MatchReaction(parsed_reaction)
        best_reaction = reaction_matches.GetBestMatch()
        
        if not best_reaction:
            return render_to_response('search_error_page.html', template_data)
        
        logging.debug('Generating a reaction from the matched KEGG IDs')
        rxn = reaction.Reaction.FromIds(best_reaction,
                                        cond=conditions.StandardConditions(),
                                        pH=ph, pMg=pmg,
                                        ionic_strength=ionic_strength)
        
        template_data = rxn.GetTemplateData(query)
        return render_to_response('reaction_page.html', template_data)

    else:
        # Otherwise we try to parse it as a single compound.
        logging.debug('Parsing the query as a single compound')
        results = matcher.Match(query)
        compound_matches = [m for m in results if m.IsCompound()]
        enzyme_matches = [m for m in results if m.IsEnzyme()]
        template_data['compound_results'] = compound_matches
        template_data['enzyme_results'] = enzyme_matches
        
        
        enzymes_first = False
        if results and results[0].IsEnzyme():
            enzymes_first = True
        template_data['enzymes_first'] = enzymes_first
            
        return render_to_response('search_results.html',
                                  template_data)

    raise Http404
