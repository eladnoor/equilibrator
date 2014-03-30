import logging

from django.http import HttpResponseBadRequest
from django.shortcuts import render_to_response
from gibbs import reaction
from gibbs import reaction_form


def HalfReactionPage(request):    
    """Renders a page for a particular reaction."""
    form = reaction_form.ReactionForm(request.GET)
    if not form.is_valid():
        logging.error(form.errors)
        return HttpResponseBadRequest('Invalid reaction form.')

    rxn = reaction.Reaction.FromForm(form)
    query = form.cleaned_query
    if form.cleaned_reactionId:
        query = rxn.GetQueryString()
    if form.cleaned_submit == 'Reverse':
        rxn.SwapSides()
        query = rxn.GetQueryString()
    
    rxn.StandardizeHalfReaction()
        
    # Render the template.
    template_data = rxn.GetTemplateData(query)
    template_data['ph'] = rxn.ph
    template_data['pmg'] = rxn.pmg
    template_data['ionic_strength'] = rxn.i_s
    return render_to_response('half_reaction_page.html', template_data)