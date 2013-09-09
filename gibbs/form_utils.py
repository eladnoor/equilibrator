from django import forms
import logging


class BaseForm(forms.Form):
    """Base form class with nice helpers."""
    
    def _GetWithDefault(self, key, default):
        if (key not in self.cleaned_data or
            self.cleaned_data[key] == None):
            return default
        return self.cleaned_data[key]
    

class ListFormField(forms.MultipleChoiceField):
    """A form field for a list of values that are unchecked.
    
    The Django MultipleChoiceField does *almost* what we want, except
    it validates that each choice is in a supplied list of choices, 
    even when that list is empty. We simply override the validation.
    """
    
    def valid_value(self, value):
        return True