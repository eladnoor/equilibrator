import logging

from django.http import HttpResponseBadRequest
from django.shortcuts import render_to_response
from gibbs import reaction
from gibbs import reaction_form


_REACTION_TEMPLATES_BY_SUBMIT = {'Update': 'reaction_page.html',
                                 'Save': 'print_reaction.html',
                                 'Reverse': 'reaction_page.html'}


def ReactionPage(request):    
    """Renders a page for a particular reaction."""
    form = reaction_form.ReactionForm(request.GET)
    if not form.is_valid():
        logging.error(form.errors)
        return HttpResponseBadRequest('Invalid reaction form.')
    
    # Figure out which template to render (based on which submit button was
    # pressed).
    rxn = reaction.Reaction.FromForm(form)
    if form.cleaned_submit == 'Reverse':
        rxn.SwapSides()
    query = rxn.GetQueryString()
    
    # Render the template.
    template_name = _REACTION_TEMPLATES_BY_SUBMIT[form.cleaned_submit]
    return render_to_response(template_name, rxn.GetTemplateData(query))