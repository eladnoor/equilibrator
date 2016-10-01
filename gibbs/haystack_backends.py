from django.conf import settings
from xapian_backend import XapianEngine, XapianSearchBackend
from xapian_backend import NGRAM_MIN_LENGTH, NGRAM_MAX_LENGTH



class CustomXapianBackend(XapianSearchBackend):

    SETTINGS = {'min_ngram_length': NGRAM_MIN_LENGTH,
                'max_ngram_length': NGRAM_MAX_LENGTH}

    def __init__(self, connection_alias, **connection_options):
        super(CustomXapianBackend, self).__init__(
            connection_alias, **connection_options)
        user_settings = getattr(settings, 'XAPIAN_SETTINGS', {})
        if user_settings:
            self.SETTINGS.update(user_settings)

        self.min_ngram_length = self.SETTINGS['min_ngram_length']
        self.max_ngram_length = self.SETTINGS['max_ngram_length']
        

    def _get_ngram_lengths(value):
        """Use the configurable n-gram length parameters.
        
        Allows for some minor memory savings.
        """
        min_l = self.min_ngram_length
        max_l = self.max_ngram_length + 1
        values = value.split()
        for item in values:
            for ngram_length in six.moves.range(min_l, max_l):
                yield item, ngram_length
                

class CustomXapianEngine(XapianEngine):
    backend = CustomXapianBackend
    