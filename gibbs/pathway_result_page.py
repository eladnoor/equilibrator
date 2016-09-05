#!/usr/bin/python

import logging
import io

from django.http import HttpResponseBadRequest
from django.http import HttpResponse
from django.shortcuts import render_to_response
from django.views.decorators.csrf import csrf_exempt
from gibbs.conditions import AqueousParams
from gibbs.forms import BuildPathwayModelForm, AnalyzePathwayModelForm
from os import path
from pathways import ParsedPathway, PathwayParseError
from pathways.bounds import Bounds
from pathways.concs import ConcentrationConverter, NoSuchUnits
from util.SBtab import SBtabTools

RELPATH = path.dirname(path.realpath(__file__))
COFACTORS_FNAME = path.join(RELPATH, '../pathways/data/cofactors.csv')


def make_aq_params(form):
    return AqueousParams(
            pH=form.cleaned_data['pH'],
            ionic_strength=form.cleaned_data['ionic_strength'])


def make_bounds(form):
    max_c = form.cleaned_data.get('max_c')
    min_c = form.cleaned_data.get('min_c')
    bounds_units = form.cleaned_data.get('conc_units')

    max_c = ConcentrationConverter.to_molar_string(max_c, bounds_units)
    min_c = ConcentrationConverter.to_molar_string(min_c, bounds_units)

    bounds = Bounds.from_csv_filename(
        COFACTORS_FNAME, default_lb=min_c, default_ub=max_c)
    return bounds


def BuildPathwayModel(request):
    """Renders a page for a particular reaction."""
    form = BuildPathwayModelForm(request.POST, request.FILES)
    if not form.is_valid():
        logging.error(form.errors)
        return HttpResponseBadRequest('Invalid pathway form.')

    try:
        bounds = make_bounds(form)
        aq_params = make_aq_params(form)
    except Exception as e:
        logging.error(e)
        return HttpResponseBadRequest(e)

    # TODO handle custom concentration bounds
    try:
        f = request.FILES['pathway_file']
        fname_base, ext = path.splitext(f.name)
        output_fname = '%s_pH%.2f_I%.2f.tsv' % (
            fname_base, aq_params.pH, aq_params.ionic_strength)
        logging.info(output_fname)

        f_data = unicode(f.read())
        sio = io.StringIO(f_data, newline=None)  # universal newline mode
        pp = ParsedPathway.from_csv_file(sio)
    except PathwayParseError as ppe:
        logging.error(ppe)
        return HttpResponseBadRequest(ppe)

    response = HttpResponse(content_type='text/tab-separated-values')
    response['Content-Disposition'] = 'attachment; filename="%s"' % output_fname
    response.write(pp.to_full_sbtab())

    return response


def read_sbtabs(f):
    """Return reactions, fluxes, keqs, bounds."""
    sbtabs = SBtabTools.openMultipleSBtabFromFile(f)
    tdict = dict([(t.getTableInformation()[1].upper(), t) for t in sbtabs])
    expected_tnames = ['REACTION', 'RELATIVEFLUX', 'REACTIONCONSTANT',
                       'CONCENTRATIONCONSTRAINT']
    assert set(expected_tnames).issubset(tdict.keys())

    return [tdict[n] for n in expected_tnames]


def PathwayResultPage(request):
    """Renders a page for a particular reaction."""
    form = AnalyzePathwayModelForm(request.POST, request.FILES)
    if not form.is_valid():
        logging.error(form.errors)
        return HttpResponseBadRequest('Invalid pathway form.')

    try:
        f_data = unicode(request.FILES['pathway_file'].read())
        sio = io.StringIO(f_data, newline=None)  # universal newline mode
        reactions, fluxes, keqs, bounds = read_sbtabs(sio)
        pp = ParsedPathway.from_full_sbtab(
            reactions, fluxes, bounds, keqs)
        logging.info('Parsed pathway.')
    except PathwayParseError as ppe:
        logging.error(ppe)
        return HttpResponseBadRequest(ppe.message)

    try:
        # calculate the MDF with the specified bounds. Render template.
        mdf_result = pp.calc_mdf()
        template_data = {'pathway': pp,
                         'mdf_result': mdf_result}
        logging.info('Calculated MDF %s', mdf_result.mdf)
    except Exception as e:
        logging.error(e)
        return HttpResponseBadRequest(e.message)
    return render_to_response('pathway_result_page.html', template_data)
