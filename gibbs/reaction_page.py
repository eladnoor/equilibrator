import logging

from django.http import HttpResponseBadRequest
from django.shortcuts import render_to_response
from gibbs import reaction
from gibbs import reaction_form


_REACTION_TEMPLATES_BY_SUBMIT = {'Update': 'reaction_page.html',
                                 'Save': 'print_reaction.html'}


def ReactionPage(request):    
    """Renders a page for a particular reaction."""
    form = reaction_form.ReactionForm(request.GET)
    if not form.is_valid():
        logging.error(form.errors)
        return HttpResponseBadRequest('Invalid reaction form.')
    
    # Figure out which template to render (based on which submit button they
    # pressed).
    template_name = _REACTION_TEMPLATES_BY_SUBMIT.get(form.cleaned_submit,
                                                      'reaction_page.html')

    rxn = reaction.Reaction.FromForm(form)
    if form.cleaned_reactionId:
        query = rxn.GetQueryString()
    elif form.cleaned_submit == 'Reverse':
        rxn.SwapSides()
        query = rxn.GetQueryString()
    elif form.cleaned_replace_co2:
        rxn.TryReplaceCO2()
        query = rxn.GetQueryString()
    else:
        query = form.cleaned_query
    
    # Render the template.
    template_data = rxn.GetTemplateData(query)
    return render_to_response(template_name, template_data)