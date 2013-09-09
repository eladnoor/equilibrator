from django.core.management import setup_environ
import settings

def SetupDjango():
    setup_environ(settings)
    
def IsSuperUser(request):
    """Returns true if the current user is a super user."""
    return request.user and request.user.is_superuser
    