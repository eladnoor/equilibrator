from django import forms

class BuildPathwayModelForm(forms.Form):
    pathway_file = forms.FileField(required=True)
    min_c = forms.FloatField(required=False)
    max_c = forms.FloatField(required=False)
    pH = forms.FloatField(required=False)
    ionic_strength = forms.FloatField(required=False)
    conc_units = forms.CharField(required=True)

class AnalyzePathwayModelForm(forms.Form):
    ANALYSIS_CHOICES = (('MDF', 'Max-min Driving Force'),
                        ('ECM', 'Enzyme Cost Minimization'))

    pathway_file = forms.FileField(required=True)
    optimization_method = forms.ChoiceField(choices=ANALYSIS_CHOICES,
                                            required=True,
                                            label='optimization_method',
                                            widget=forms.RadioSelect())
    
#    def __init__(self, custom_choices=None, *args, **kwargs):
#        super(AnalyzePathwayModelForm, self).__init__(*args, **kwargs)
#        if custom_choices:
#            self.fields['field'].choices = custom_choices