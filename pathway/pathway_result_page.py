#!/usr/bin/python

from gibbs.conditions import AqueousParams
from .bounds import Bounds
from .concs import ConcentrationConverter
from util.SBtab import SBtabTools
import os

from equilibrator.settings import BASE_DIR
COFACTORS_FNAME = os.path.join(BASE_DIR, 'pathway/data/cofactors.csv')


def make_aq_params(form):
    return AqueousParams(
            pH=form.cleaned_data['pH'],
            ionic_strength=form.cleaned_data['ionic_strength'])


def make_bounds(form):
    max_c = form.cleaned_data.get('max_c')
    min_c = form.cleaned_data.get('min_c')
    bounds_units = form.cleaned_data.get('conc_units')

    min_c = ConcentrationConverter.to_molar_string(min_c, bounds_units)
    max_c = ConcentrationConverter.to_molar_string(max_c, bounds_units)

    bounds = Bounds.from_csv_filename(
        COFACTORS_FNAME, default_lb=min_c, default_ub=max_c)
    return bounds


def read_sbtabs_mdf(f):
    """Return reactions, fluxes, keqs, bounds."""
    sbtabs = SBtabTools.openMultipleSBtabFromFile(f)
    tdict = dict([(t.getTableInformation()[1].upper(), t) for t in sbtabs])
    expected_tnames = ['REACTION', 'RELATIVEFLUX', 'REACTIONCONSTANT',
                       'CONCENTRATIONCONSTRAINT']
    assert set(expected_tnames).issubset(tdict.keys())

    return [tdict[n] for n in expected_tnames]

def read_sbtabs_ecm(f):
    """Return reactions, fluxes, keqs, bounds."""
    sbtabs = SBtabTools.openMultipleSBtabFromFile(f)
    tdict = dict([(t.getTableInformation()[1].upper(), t) for t in sbtabs])
    expected_tnames = ['REACTION', 'RELATIVEFLUX', 'REACTIONCONSTANT',
                       'CONCENTRATIONCONSTRAINT']
    assert set(expected_tnames).issubset(tdict.keys())

    return [tdict[n] for n in expected_tnames]
