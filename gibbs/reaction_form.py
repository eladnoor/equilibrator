from django import forms
from gibbs import form_utils
from gibbs import constants
from gibbs import conditions

class ReactionForm(form_utils.BaseForm):
    
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
    
    reactionId = forms.CharField(required=False)
    reactantsId = form_utils.ListFormField(required=False)
    reactantsCoeff = form_utils.ListFormField(required=False)
    reactantsName = form_utils.ListFormField(required=False)
    reactantsPhase = form_utils.ListFormField(required=False)
    reactantsConcentration = form_utils.ListFormField(required=False)

    ph = forms.FloatField(required=False)
    pmg = forms.FloatField(required=False)
    ionic_strength = forms.FloatField(required=False)
    electronReductionPotential = forms.FloatField(required=False)
    conditions = forms.ChoiceField(required=False,
                                   choices=constants.CONDITION_CHOICES)
    
    query = forms.CharField(max_length=2048, required=False)
    balance_w_water = forms.BooleanField(required=False)
    balance_electrons = forms.BooleanField(required=False)
    replace_co2 = forms.BooleanField(required=False)
    submit = forms.ChoiceField(required=False,
                               choices=[('Update', 'update'),
                                        ('Save', 'save'),
                                        ('Reverse', 'reverse')])
    
    # Convenience accessors for clean data with defaults.
    cleaned_reactionId = property(lambda self: self.cleaned_data['reactionId'])
    cleaned_reactantsId = property(lambda self: self.cleaned_data['reactantsId'])
    cleaned_reactantsCoeff = property(lambda self: [float(c) for c in self.cleaned_data['reactantsCoeff']])
    cleaned_reactantsName = property(lambda self: self.cleaned_data['reactantsName'])
    cleaned_reactantsPhase = property(lambda self: self.cleaned_data['reactantsPhase'])
    cleaned_reactantsConcentration = property(GetReactantConcentrations)
    cleaned_ph = property(lambda self: self._GetWithDefault('ph', constants.DEFAULT_PH))
    cleaned_pmg = property(lambda self: self._GetWithDefault('pmg', constants.DEFAULT_PMG))
    cleaned_ionic_strength = property(lambda self: self._GetWithDefault('ionic_strength',
                                                    constants.DEFAULT_IONIC_STRENGTH))
    cleaned_e_reduction_potential = property(lambda self: self._GetWithDefault('electronReductionPotential',
                                                    constants.DEFAULT_ELECTRON_REDUCTION_POTENTIAL))
    cleaned_conditions = property(lambda self: self.cleaned_data['conditions'])
    cleaned_query = property(lambda self: self.cleaned_data['query'])
    cleaned_balance_w_water = property(lambda self: self._GetWithDefault('balance_w_water', False))
    cleaned_balance_electrons = property(lambda self: self._GetWithDefault('balance_electrons', False))
    cleaned_replace_co2 = property(lambda self: self._GetWithDefault('replace_co2', False))
    cleaned_submit = property(lambda self: self._GetWithDefault('submit', 'Update'))
