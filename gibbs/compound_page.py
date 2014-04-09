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
    
    rxn = reaction.Reaction.FromForm(form)

    if len(rxn.reactants) != 1:
        logging.error('There must be only 1 reactant in a "compound" page')
        raise Http404

    rxn.aq_params = conditions.AqueousParams.FromForm(form, request.COOKIES) 
    logging.info('Aqueous parameters: ' + str(rxn.aq_params))
    
    compound = rxn.reactants[0].compound
    compound.StashTransformedSpeciesEnergies(pH=rxn.pH,
                                             pMg=rxn.pMg,
                                             ionic_strength=rxn.ionic_strength)
    
    template_data = rxn.GetTemplateData(None)
    template_data.update({'compound': rxn.reactants[0].compound})
    response = render_to_response('compound_page.html', template_data)
    rxn.aq_params.SetCookies(response)
    return response