from haystack import indexes
from django.apps import apps
from gibbs.models import CommonName


class CommonNameIndex(indexes.SearchIndex, indexes.Indexable):
    """
        Index for searching CommonNames.
    """
    text = indexes.CharField(document=True, use_template=True)

    # We add this for autocomplete.
    title_autocomplete = indexes.NgramField(model_attr='name')

    def get_model(self):
        return CommonName

    def index_queryset(self, using=None):
        """Used when the entire index for model is updated."""
        return apps.get_model('gibbs.CommonName').objects.all()
