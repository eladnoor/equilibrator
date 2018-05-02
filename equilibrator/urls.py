"""equilibrator URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/2.0/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.urls import path, include
from django.views.generic import TemplateView
from gibbs import views

urlpatterns = [
    path('admin/', admin.site.urls),
#    path('gibbs/', include('gibbs.urls')),
    path('', TemplateView.as_view(template_name='main.html'), name='index'),
    path(r'compound_image', views.CompoundImage),
    path(r'compound', views.CompoundPage),
    path(r'download', views.DownloadPage),
    path(r'enzyme', views.EnzymePage),
    path(r'reaction', views.ReactionPage),
    path(r'graph_reaction', views.ReactionGraph),
    path(r'data_refs', views.RefsPage),
    path(r'search', views.ResultsPage),
    path(r'suggest', views.SuggestJson),
    path('pathway/', include('pathway.urls')),
    path(r'robots\.txt', TemplateView.as_view(template_name='robots.txt')),
] 
