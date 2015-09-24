from gibbs.models import CommonName
from haystack import indexes


class CommonNameIndex(indexes.SearchIndex, indexes.Indexable):
    """Index for searching CommonNames."""
    text = indexes.CharField(document=True, use_template=True)
    
    # We add this for autocomplete.
    name_auto = indexes.NgramField(model_attr='name')

    def get_model(self):
        return CommonName

    def index_queryset(self, using=None):
        """Used when the entire index for model is updated."""
        return CommonName.objects.all()