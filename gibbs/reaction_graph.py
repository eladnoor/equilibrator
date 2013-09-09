import json
import logging

from gibbs import reaction
from gibbs import reaction_graph_form
from django.shortcuts import render_to_response


def ReactionGraph(request):
    """Renders the graph page."""
    form = reaction_graph_form.ReactionGraphForm(request.GET)
    if not form.is_valid():
        logging.error(form.errors)
        raise Http404

    mode = 'varyPh'
    if form.cleaned_vary_is:
        mode = 'varyIs'

    rxn = reaction.Reaction.FromForm(form)
    compound_list = [c.compound.SpeciesJson()['species'] for c in rxn.filtered_substrates] + \
                    [c.compound.SpeciesJson()['species'] for c in rxn.filtered_products]
    coeff_list = [-c.coeff for c in rxn.filtered_substrates] + \
                 [c.coeff for c in rxn.filtered_products]
    concentration_list = [c.concentration for c in rxn.filtered_substrates] + \
                         [c.concentration for c in rxn.filtered_products]
    template_data = {"compound_data": json.dumps(compound_list),
                     "coeff_data": json.dumps(coeff_list),
                     "concentration_list": json.dumps(concentration_list),
                     "mode": mode,
                     "reaction": rxn}
    
    return render_to_response('reaction_graph.html', template_data)