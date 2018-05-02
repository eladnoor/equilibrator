from django.urls import path
from django.contrib import admin
from . import views
admin.autodiscover()

urlpatterns = [
    path(r'', views.DefinePathwayPage, name='index'),
    path(r'build_model', views.BuildPathwayModel),
    path(r'results', views.PathwayResultPage),
]
