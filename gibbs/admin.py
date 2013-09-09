import inspect

from gibbs import models
from django.contrib import admin

admin.site.register(models.CommonName)
admin.site.register(models.Compound)
admin.site.register(models.Specie)
admin.site.register(models.ValueSource)
admin.site.register(models.StoredReaction)
admin.site.register(models.Enzyme)
