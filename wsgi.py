import os
import sys

os.environ['DJANGO_SETTINGS_MODULE'] = 'equilibrator.settings'

from django.core.wsgi import get_wsgi_application

application = get_wsgi_application()
