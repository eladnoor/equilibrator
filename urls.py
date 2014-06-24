from django.conf.urls import patterns, include
from django.conf import settings

# Uncomment the next two lines to enable the admin:
from django.contrib import admin
admin.autodiscover()

urlpatterns = patterns('',
    (r'^media/(?P<path>.*)$', 'django.views.static.serve', {'document_root': settings.MEDIA_ROOT }),
    (r'^$', 'equilibrator.gibbs.main_page.MainPage'),
    (r'^about', 'equilibrator.gibbs.info_pages.AboutPage'),
    (r'^faq', 'equilibrator.gibbs.info_pages.FAQPage'),
    (r'^whats_new', 'equilibrator.gibbs.info_pages.WhatsNewPage'),
    (r'^redox_review', 'equilibrator.gibbs.info_pages.RedoxReview'),
    (r'^cite', 'equilibrator.gibbs.info_pages.CitePage'),
    (r'^classic_reactions', 'equilibrator.gibbs.info_pages.ClassicReactions'),
    (r'^compound_data', 'equilibrator.gibbs.compound_api.CompoundAPI'),
    (r'^compound_image', 'equilibrator.gibbs.compound_image.CompoundImage'),
    (r'^compound', 'equilibrator.gibbs.compound_page.CompoundPage'),
    (r'^download', 'equilibrator.gibbs.info_pages.DownloadPage'),
    (r'^enzyme', 'equilibrator.gibbs.enzyme_page.EnzymePage'),
    (r'^reaction', 'equilibrator.gibbs.reaction_page.ReactionPage'),
    (r'^graph_reaction', 'equilibrator.gibbs.reaction_graph.ReactionGraph'),
    (r'^data_refs', 'equilibrator.gibbs.data_refs_page.RefsPage'),
    (r'^search', 'equilibrator.gibbs.search_results_page.ResultsPage'),
    (r'^suggest', 'equilibrator.gibbs.suggest.SuggestJson'),
    (r'^robots\.txt', 'equilibrator.gibbs.info_pages.Robots'),
    # Example:
    # (r'^equilibrator/', include('equilibrator.foo.urls')),

    # Uncomment the admin/doc line below and add 'django.contrib.admindocs' 
    # to INSTALLED_APPS to enable admin documentation:
    # (r'^admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
    (r'^admin/', include(admin.site.urls)),
)
