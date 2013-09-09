from django import forms
from gibbs import constants
from gibbs import form_utils


class CompoundForm(form_utils.BaseForm):
    compoundId = forms.CharField(max_length=50)
    ph = forms.FloatField(required=False)
    pmg = forms.FloatField(required=False)
    ionic_strength = forms.FloatField(required=False)
    
    # Convenience accessors for clean data with defaults.
    cleaned_compoundId = property(lambda self: self.cleaned_data['compoundId'])
    cleaned_ph = property(lambda self: self._GetWithDefault('ph',
                                                            constants.DEFAULT_PH))
    cleaned_pmg = property(lambda self: self._GetWithDefault('pmg', constants.DEFAULT_PMG))
    cleaned_ionic_strength = property(lambda self: self._GetWithDefault('ionic_strength',
                                                                        constants.DEFAULT_IONIC_STRENGTH))