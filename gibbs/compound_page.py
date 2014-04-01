import logging

from django.http import Http404
from django.shortcuts import render_to_response
from gibbs import compound_form
from gibbs import reaction

def CompoundPage(request):
    """Renders a page for a particular compound."""
    form = compound_form.CompoundForm(request.GET)
    if not form.is_valid():
        logging.error(form.errors)
        raise Http404
    
    rxn = reaction.Reaction.FromForm(form)    
    compound = rxn.reactants[0].compound
    compound.StashTransformedSpeciesEnergies(pH=form.cleaned_ph,
                                             pMg=form.cleaned_pmg,
                                             ionic_strength=form.cleaned_ionic_strength)
    
    template_data = rxn.GetTemplateData(None)
    template_data.update({'compound': rxn.reactants[0].compound})
    return render_to_response('compound_page.html', template_data)