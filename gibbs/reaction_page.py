import logging

from django.http import HttpResponseBadRequest, Http404
from django.shortcuts import render_to_response
from gibbs import reaction
from gibbs.forms import ReactionForm
from gibbs import conditions

_REACTION_TEMPLATES_BY_SUBMIT = {'': 'reaction_page.html',
                                 'Update': 'reaction_page.html',
                                 'Save': 'print_reaction.html',
                                 'Reverse': 'reaction_page.html',
                                 'Reset': 'reaction_page.html'}


def ReactionPage(request):    
    """Renders a page for a particular reaction."""
    form = ReactionForm(request.GET)
    if not form.is_valid():
        logging.error(form.errors)
        return HttpResponseBadRequest('Invalid reaction form.')
    
    rxn = reaction.Reaction.FromForm(form)

    rxn.aq_params = conditions.AqueousParams.FromForm(form, request.COOKIES) 
    logging.info('Aqueous parameters: ' + str(rxn.aq_params))

    if form.cleaned_submit == 'Reverse':
        rxn.SwapSides()
    elif form.cleaned_submit == 'Reset':
        rxn.aq_params = conditions.AqueousParams()
        rxn.ResetConcentrations()
    query = rxn.GetQueryString()
    
    # Render the template.
    if form.cleaned_submit not in _REACTION_TEMPLATES_BY_SUBMIT:
        logging.error('Unknown submit term for reaction page: ' + form.cleaned_submit)
        raise Http404
    template_name = _REACTION_TEMPLATES_BY_SUBMIT[form.cleaned_submit]
    response = render_to_response(template_name, rxn.GetTemplateData(query))    
    rxn.aq_params.SetCookies(response)
    return response
