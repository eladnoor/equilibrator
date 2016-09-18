import os

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "settings")
from django.conf import settings

import logging

logger = logging.getLogger('')
logger.setLevel(logging.DEBUG)


def IsSuperUser(request):
    """Returns true if the current user is a super user."""
    return request.user and request.user.is_superuser
