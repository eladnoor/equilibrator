import os

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "settings")
from django.conf import settings
    
def IsSuperUser(request):
    """Returns true if the current user is a super user."""
    return request.user and request.user.is_superuser
    