import views
from django.conf import settings
from django.conf.urls import include, url
from django.conf.urls.static import static
from django.contrib import admin
admin.autodiscover()

urlpatterns = [
    url(r'^$', views.MainPage),
    url(r'^about', views.AboutPage),
    url(r'^faq', views.FAQPage),
    url(r'^whats_new', views.WhatsNewPage),
    url(r'^cite', views.CitePage),
    url(r'^classic_reactions', views.ClassicReactions),
    url(r'^compound_image', views.CompoundImage),
    url(r'^compound', views.CompoundPage),
    url(r'^download', views.DownloadPage),
    url(r'^enzyme', views.EnzymePage),
    url(r'^reaction', views.ReactionPage),
    url(r'^graph_reaction', views.ReactionGraph),
    url(r'^data_refs', views.RefsPage),
    url(r'^search', views.ResultsPage),
    url(r'^suggest', views.SuggestJson),
    url(r'^pathway$', views.DefinePathwayPage),
    url(r'^pathway/build_model', views.BuildPathwayModel),
    url(r'^pathway/results', views.PathwayResultPage),
    url(r'^robots\.txt', views.Robots),
    # Example:
    # (r'^equilibrator/', include('equilibrator.foo.urls')),

    # Uncomment the admin/doc line below and add 'django.contrib.admindocs'
    # to INSTALLED_APPS to enable admin documentation:
    # (r'^admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
    url(r'^admin/', include(admin.site.urls)),
]
