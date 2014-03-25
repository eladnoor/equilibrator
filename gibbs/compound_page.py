import logging

from django.http import Http404
from django.shortcuts import render_to_response
from gibbs import compound_form
from gibbs import models
from util import django_utils


def CompoundPage(request):
    """Renders a page for a particular compound."""
    form = compound_form.CompoundForm(request.GET)
    if not form.is_valid():
        logging.error(form.errors)
        raise Http404
    
    # Compute the delta G estimate.
    kegg_id = form.cleaned_compoundId
    compound = models.Compound.objects.get(kegg_id=kegg_id)
    compound.StashTransformedSpeciesEnergies(form.cleaned_ph,
                                             form.cleaned_pmg,
                                             form.cleaned_ionic_strength)
    delta_g0_prime_estimate = compound.DeltaG0Prime(
        pH=form.cleaned_ph,
        pMg=form.cleaned_pmg,
        ionic_strength=form.cleaned_ionic_strength)
    
    template_data = {'is_superuser': django_utils.IsSuperUser(request),
                     'compound': compound,
                     'ph': form.cleaned_ph,
                     'pmg': form.cleaned_pmg,
                     'ionic_strength': form.cleaned_ionic_strength,
                     'delta_g0_prime_estimate': delta_g0_prime_estimate,
                     'no_dg_explanation': compound.no_dg_explanation,
                     'kegg_link': compound.GetKeggLink()}
    return render_to_response('compound_page.html', template_data)