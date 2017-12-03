from django import forms

class BuildPathwayModelForm(forms.Form):
    pathway_file = forms.FileField(required=True)
    min_c = forms.FloatField(required=False)
    max_c = forms.FloatField(required=False)
    pH = forms.FloatField(required=False)
    ionic_strength = forms.FloatField(required=False)
    conc_units = forms.CharField(required=True)

class AnalyzePathwayModelForm(forms.Form):
    pathway_file = forms.FileField(required=True)
    pH = forms.FloatField(required=False)
    ionic_strength = forms.FloatField(required=False)
