# Django settings for DjangoTest project.

import os
import sys

import matplotlib

PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.normcase(PROJECT_ROOT))
sys.path.insert(0, os.path.join(PROJECT_ROOT, ".."))

# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

DEBUG = True

ALLOWED_HOSTS = [
    '*',  # Allow domain and subdomains
]

ADMINS = (
    ('Avi Flamholz', 'flamholz@gmail.com'),
    ('Elad Noor', 'elad.noor@gmail.com'),
)

MANAGERS = ADMINS

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.mysql',
        'NAME': 'milolab_eqbtr',
        'USER': 'milolab_eqbtr',
        'PASSWORD': 'password',
        'HOST': '',
        'PORT': '',
    }
}

# Local time zone for this installation. Choices can be found here:
# http://en.wikipedia.org/wiki/List_of_tz_zones_by_name
# although not all choices may be available on all operating systems.
# On Unix systems, a value of None will cause Django to use the same
# timezone as the operating system.
# If running in a Windows environment this must be set to the same as your
# system time zone.
TIME_ZONE = 'America/Chicago'

# Language code for this installation. All choices can be found here:
# http://www.i18nguy.com/unicode/language-identifiers.html
LANGUAGE_CODE = 'en-us'

SITE_ID = 1

# If you set this to False, Django will make some optimizations so as not
# to load the internationalization machinery.
USE_I18N = True

# If you set this to False, Django will not format dates, numbers and
# calendars according to the current locale
USE_L10N = True

# Absolute path to the directory that holds media.
# Example: "/home/media/media.lawrence.com/"
MEDIA_ROOT = os.path.join(PROJECT_ROOT, 'media')

# URL that handles the media served from MEDIA_ROOT. Make sure to use a
# trailing slash if there is a path component (optional in other cases).
# Examples: "http://media.lawrence.com", "http://example.com/media/"
MEDIA_URL = 'http://127.0.0.1:8000/media/'

# URL prefix for admin media -- CSS, JavaScript and images. Make sure to use a
# trailing slash.
# Examples: "http://foo.com/media/", "/media/".
ADMIN_MEDIA_PREFIX = '/admin_media/'

# for serving static files
STATIC_URL = '/static/'
STATICFILES_DIRS = [os.path.join(PROJECT_ROOT, "static")]

# Make this unique, and don't share it with anybody.
SECRET_KEY = 'b88!&88r-69=r*%q8cgnj&9dfm!^1u!ij3+jnkoebh4vrm41we'

MIDDLEWARE_CLASSES = (
    'django.middleware.common.CommonMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
)

ROOT_URLCONF = 'equilibrator.urls'

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [os.path.join(PROJECT_ROOT, 'templates')],
        'APP_DIRS': True,
        'OPTIONS': {
            # the 'string_if_invalid' should not be left uncommented
            # except when debugging a spcific template issue
            'string_if_invalid': 'INVALID %s',
            'context_processors': [
                'django.template.context_processors.debug',
                'django.template.context_processors.request',
                'django.contrib.auth.context_processors.auth',
                'django.contrib.messages.context_processors.messages',
            ],
        },
    },
]

INSTALLED_APPS = (
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.sites',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'haystack',
    # Uncomment the next line to enable the admin:
    'django.contrib.admin',
    'django.contrib.humanize',
    'django_extensions',
    'gibbs',
)

# Haystack related settings - for search/autocomplete.
HAYSTACK_CONNECTIONS = {
    'default': {
                'ENGINE': 'haystack.backends.solr_backend.SolrEngine',
                'URL': 'http://127.0.0.1:8080/solr',
                },
}


# Custom Xapian settings
XAPIAN_SETTINGS = {
    'min_ngram_length': 2,
    'max_ngram_length': 8,
}

SESSION_ENGINE = 'django.contrib.sessions.backends.signed_cookies'

LOGGING = {
    'version': 1,
    'disable_existing_loggers': True,
    'formatters': {
        'verbose': {
            'format': "[%(asctime)s %(filename)s:%(lineno)s] %(levelname)s %(message)s",
            'datefmt': "%d/%b/%Y %H:%M:%S"
        },
        'simple': {
            'format': '%(levelname)s %(message)s'
        },
    },
    'handlers': {
        'file': {
            'level': 'INFO',
            'class': 'logging.FileHandler',
            'filename': os.path.join(PROJECT_ROOT, 'gibbs.log'),
            'formatter': 'verbose'
        },
        'console': {
            'level': 'INFO',
            'class': 'logging.StreamHandler',
            'formatter': 'verbose'
        },
    },
    'loggers': {
        '': {
            'handlers': ['console'],
            'level': 'INFO',
        },
    }
}

# DJANGO-PROFILER 2.0
PROFILING_LOGGER_NAME = 'profiler.log'
PROFILING_SQL_QUERIES = False
