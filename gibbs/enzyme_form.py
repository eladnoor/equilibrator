from django import forms
from gibbs import form_utils


class EnzymeForm(form_utils.BaseForm):
    ec = forms.CharField(max_length=50)
    
    # Convenience accessors for clean data with defaults.
    cleaned_ec = property(lambda self: self.cleaned_data['ec'])