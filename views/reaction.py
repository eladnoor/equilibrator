import logging
import json
from django.http import HttpResponseBadRequest
from django.shortcuts import render
from django.http import Http404
from django.apps import apps
from ..gibbs import conditions
from ..gibbs.forms import ReactionForm, ReactionGraphForm

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
    aq_params = conditions.AqueousParams.FromForm(form, request.COOKIES)
    rxn = apps.get_model('gibbs.Reaction').FromForm(form, aq_params)

    if form.cleaned_submit == 'Reverse':
        rxn.SwapSides()
    elif form.cleaned_submit == 'Reset':
        rxn.aq_params = conditions.AqueousParams()
        rxn.ResetConcentrations()
    query = rxn.GetQueryString()
    # Render the template.
    if form.cleaned_submit not in _REACTION_TEMPLATES_BY_SUBMIT:
        logging.error('Unknown submit term for reaction page: ' +
                      form.cleaned_submit)
        raise Http404
    template_name = _REACTION_TEMPLATES_BY_SUBMIT[form.cleaned_submit]
    response = render(request, template_name, rxn.GetTemplateData(query))
    rxn.aq_params.SetCookies(response)
    return response


def ReactionGraph(request):
    """Renders the graph page."""
    form = ReactionGraphForm(request.GET)
    if not form.is_valid():
        logging.error(form.errors)
        raise Http404

    mode = 'varyPh'
    if form.cleaned_vary_is:
        mode = 'varyIs'

    aq_params = conditions.AqueousParams.FromForm(form, request.COOKIES)
    rxn = apps.get_model('gibbs.Reaction').FromForm(form, aq_params)

    compound_list = [c.GetCompoundList() for c in rxn.reactants]
    coeff_list = [c.coeff for c in rxn.reactants]
    concentration_list = [c.phase.Value() for c in rxn.reactants]
    template_data = {'compound_data': json.dumps(compound_list),
                     'coeff_data': json.dumps(coeff_list),
                     'concentration_list': json.dumps(concentration_list),
                     'mode': mode,
                     'reaction': rxn,
                     'query': rxn.GetQueryString()}
    template_data.update(rxn.aq_params.GetTemplateData())
    return render(request, 'reaction_graph.html', template_data)
