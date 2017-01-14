import base64
import logging
import os
from django.http import HttpResponse, HttpResponseBadRequest
from django.views.decorators.csrf import csrf_exempt
from django.shortcuts import render
from django.http import Http404
from django.apps import apps
from gibbs import conditions
from gibbs.forms import CompoundForm
from settings import STATICFILES_DIRS

NO_STRUCTURE_THUMBNAIL = os.path.join(STATICFILES_DIRS[0], 'images',
                                      'structure_not_available.png')


@csrf_exempt
def CompoundImage(request):
    data = request.GET
    compound_id = data.get('compoundId', None)
    if compound_id is None:
        return HttpResponseBadRequest('No request data.')
    compounds = apps.get_model('gibbs.Compound').objects.filter(
        kegg_id=compound_id)
    if not compounds:
        return HttpResponseBadRequest('No such compound.')
    compound = compounds[0]
    if not compound.thumbnail or compound.thumbnail == 'error':
        logging.info('No structure available for %s' % compound)
        image_data = open(NO_STRUCTURE_THUMBNAIL, 'r').read()

    else:
        image_data = base64.decodestring(compound.thumbnail)
    return HttpResponse(image_data, content_type='image/png')


def CompoundPage(request):
    """Renders a page for a particular compound."""
    form = CompoundForm(request.GET)
    if not form.is_valid():
        logging.error(form.errors)
        raise Http404
    aq_params = conditions.AqueousParams.FromForm(form, request.COOKIES)
    rxn = apps.get_model('gibbs.Reaction').FromForm(form, aq_params)
    if len(rxn.reactants) != 1:
        logging.error('There must be only 1 reactant in a "compound" page')
        raise Http404

    compound = rxn.reactants[0].compound
    logging.info('compoundID = ' + form.cleaned_compoundId)
    logging.info('rxn = ' + rxn.GetQueryString())
    logging.info('Submit = ' + form.cleaned_submit)

    if form.cleaned_submit == 'Reset':
        logging.info('resetting conditions')
        rxn.aq_params = conditions.AqueousParams()
        rxn.ResetConcentrations()

    compound.StashTransformedSpeciesEnergies(rxn.aq_params)

    query = compound.FirstName()
    template_data = rxn.GetTemplateData(query)
    template_data.update({'reaction': rxn,
                          'compound': compound,
                          'alberty_link': compound.GetNewPriorityLink(99),
                          'cc_link': compound.GetNewPriorityLink(1)})
    response = render(request, 'compound_page.html', template_data)
    rxn.aq_params.SetCookies(response)
    return response
