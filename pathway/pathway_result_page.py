#!/usr/bin/python

from gibbs.conditions import AqueousParams
from .bounds import Bounds
from .concs import ConcentrationConverter
from . import MaxMinDrivingForce, EnzymeCostMinimization
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

def from_csv(form, fp):
    bounds = make_bounds(form)
    aq_params = make_aq_params(form)
    optimization_method = form.GetOptimizationMethod()
    if optimization_method == 'MDF':
        pp = MaxMinDrivingForce.from_csv(fp, bounds=bounds,
                                         aq_params=aq_params)
    else:
        pp = EnzymeCostMinimization.from_csv(fp, bounds=bounds,
                                             aq_params=aq_params)

    fname_base, ext = os.path.splitext(fp.name)
    output_fname = '%s_pH%.2f_I%.2f_%s.tsv' % (
        fname_base, aq_params.pH,
        aq_params.ionic_strength, optimization_method)

    return pp, output_fname