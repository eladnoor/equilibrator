from django import forms
from gibbs import reaction_form

class ReactionGraphForm(reaction_form.ReactionForm):
    vary_ph = forms.BooleanField(required=False)
    vary_is = forms.BooleanField(required=False)
    vary_pmg = forms.BooleanField(required=False)
    
    # Convenience accessors for clean data with defaults.
    cleaned_vary_ph  = property(lambda self: self._GetWithDefault('vary_ph', False))
    cleaned_vary_pmg = property(lambda self: self._GetWithDefault('vary_pmg', False))
    cleaned_vary_is  = property(lambda self: self._GetWithDefault('vary_is', False))