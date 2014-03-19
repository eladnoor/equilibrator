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
    if form.cleaned_balance_w_water:
        rxn.TryBalanceWithWater()
        query = rxn.GetQueryString()
    if form.cleaned_replace_co2:
        rxn.TryReplaceCO2()
        query = rxn.GetQueryString()
    
    rxn.StandardizeHalfReaction()
        
    # Render the template.
    balance_with_water_link = rxn.GetBalanceWithWaterLink(query)
    replace_co2_link = rxn.GetReplaceCO2Link(query)
    template_data = {'reaction': rxn,
                     'query': query,
                     'ph': rxn.ph,
                     'pmg': rxn.pmg,
                     'ionic_strength': rxn.i_s,
                     'conditions': str(rxn.conditions),
                     'balance_with_water_link': balance_with_water_link,
                     'replace_co2_link': replace_co2_link}
    return render_to_response('half_reaction_page.html', template_data)