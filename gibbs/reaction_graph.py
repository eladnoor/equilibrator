import json
import logging

from gibbs import reaction
from gibbs.forms import ReactionGraphForm
from django.shortcuts import render_to_response
from django.http import Http404
from gibbs import conditions

def ReactionGraph(request):
    """Renders the graph page."""
    form = ReactionGraphForm(request.GET)
    if not form.is_valid():
        logging.error(form.errors)
        raise Http404

    mode = 'varyPh'
    if form.cleaned_vary_is:
        mode = 'varyIs'

    logging.info('reading reaction graph form')

    aq_params = conditions.AqueousParams.FromForm(form, request.COOKIES) 
    rxn = reaction.Reaction.FromForm(form, aq_params)

    logging.info([c.phase for c in rxn.reactants])
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
    return render_to_response('reaction_graph.html', template_data)