import logging

from django.http import Http404
from django.shortcuts import render_to_response
from gibbs.forms import CompoundForm
from gibbs import reaction
from gibbs import conditions

def CompoundPage(request):
    """Renders a page for a particular compound."""
    form = CompoundForm(request.GET)
    if not form.is_valid():
        logging.error(form.errors)
        raise Http404
    
    aq_params = conditions.AqueousParams.FromForm(form, request.COOKIES) 
    rxn = reaction.Reaction.FromForm(form, aq_params)

    if len(rxn.reactants) != 1:
        logging.error('There must be only 1 reactant in a "compound" page')
        raise Http404

    compound = rxn.reactants[0].compound

    logging.info('Submit = ' + form.cleaned_submit)

    if form.cleaned_submit == 'Reset':
        logging.info('resetting conditions')
        rxn.aq_params = conditions.AqueousParams()
        rxn.ResetConcentrations()

    compound.StashTransformedSpeciesEnergies(rxn.aq_params)

    query = compound.FirstName()
    template_data = rxn.GetTemplateData(query)
    template_data.update({'compound': compound,
                          'alberty_link': compound.GetNewPriorityLink(99),
                          'cc_link': compound.GetNewPriorityLink(1)})
    response = render_to_response('compound_page.html', template_data)
    rxn.aq_params.SetCookies(response)
    return response