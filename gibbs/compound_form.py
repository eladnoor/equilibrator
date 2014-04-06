from django import forms
from gibbs import constants
from gibbs import form_utils

class CompoundForm(form_utils.BaseForm):
   
    def GetReactantConcentrations(self):
        for c in self.cleaned_data['reactantsConcentration']:
            try:
                conc = float(c)
                if conc <= 0:
                    yield 1e-9
                else:
                    yield conc
            except ValueError:
                yield 1e-9

    compoundId = forms.CharField(max_length=50)
    ph = forms.FloatField(required=False)
    pmg = forms.FloatField(required=False)
    ionic_strength = forms.FloatField(required=False)
    reactantsPhase = form_utils.ListFormField(required=False)
    reactantsConcentration = form_utils.ListFormField(required=False)
    conditions = forms.ChoiceField(required=False,
                                   choices=constants.CONDITION_CHOICES)

    # Convenience accessors for clean data with defaults.
    cleaned_compoundId = property(lambda self: self.cleaned_data['compoundId'])
    cleaned_ph = property(lambda self: self._GetWithDefault('ph',
                                                            constants.DEFAULT_PH))
    cleaned_pmg = property(lambda self: self._GetWithDefault('pmg',
                                                             constants.DEFAULT_PMG))
    cleaned_ionic_strength = property(lambda self: self._GetWithDefault('ionic_strength',
                                                                        constants.DEFAULT_IONIC_STRENGTH))
    cleaned_reactantsPhase = property(lambda self: self.cleaned_data['reactantsPhase'])
    cleaned_reactantsConcentration = property(GetReactantConcentrations)
    cleaned_conditions = property(lambda self: self.cleaned_data['conditions'])

    # we need to create the following properties in order for this form
    # to impersonate a reaction_form (something we need to do for creating
    # a Reaction object using .FromForm(form))
    cleaned_reactantsCoeff = property(lambda self: [1])
    cleaned_reactantsId = property(lambda self: [self.cleaned_compoundId])
    cleaned_reactantsName = property(lambda self: [None])
    cleaned_reactantsPhase = property(lambda self: [None])
    cleaned_e_reduction_potential = property(lambda self: 0)
    cleaned_reactionId = property(lambda self: None)