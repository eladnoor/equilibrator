from django import forms
ANALYSIS_CHOICES = {'MDF': 'Max-min Driving Force',
                    'ECM': 'Enzyme Cost Minimization'}

class BuildPathwayModelForm(forms.Form):
    pathway_file = forms.FileField(required=True)
    min_c = forms.FloatField(required=False)
    max_c = forms.FloatField(required=False)
    pH = forms.FloatField(required=False)
    ionic_strength = forms.FloatField(required=False)
    conc_units = forms.CharField(required=True)
    mdf_ecm_toggle = forms.BooleanField(required=False)
    
    def GetOptimizationMethod(self):
        if self.cleaned_data['mdf_ecm_toggle']:
            optimization_method = 'ECM'
        else:
            optimization_method = 'MDF'
        return optimization_method

class AnalyzePathwayModelForm(forms.Form):
    pathway_file = forms.FileField(required=True)
    mdf_ecm_toggle = forms.BooleanField(required=False)

    def GetOptimizationMethod(self):
        if self.cleaned_data['mdf_ecm_toggle']:
            optimization_method = 'ECM'
        else:
            optimization_method = 'MDF'
        return optimization_method
    
#    def __init__(self, custom_choices=None, *args, **kwargs):
#        super(AnalyzePathwayModelForm, self).__init__(*args, **kwargs)
#        if custom_choices:
#            self.fields['field'].choices = custom_choices