import views.compound
import views.reaction
import views.views
from django.conf.urls import include, url
from django.contrib import admin
from django.views.generic import TemplateView
admin.autodiscover()

urlpatterns = [
    url(r'^$', TemplateView.as_view(template_name='main.html')),
    url(r'^about', TemplateView.as_view(template_name='about.html')),
    url(r'^faq', TemplateView.as_view(template_name='faq.html')),
    url(r'^whats_new', TemplateView.as_view(template_name='new_in_2_0.html')),
    url(r'^cite', TemplateView.as_view(template_name='cite.html')),
    url(r'^classic_reactions',
        TemplateView.as_view(template_name='classic_reactions.html')),
    url(r'^compound_image', views.compound.CompoundImage),
    url(r'^compound', views.compound.CompoundPage),
    url(r'^download', views.views.DownloadPage),
    url(r'^enzyme', views.views.EnzymePage),
    url(r'^reaction', views.reaction.ReactionPage),
    url(r'^graph_reaction', views.reaction.ReactionGraph),
    url(r'^data_refs', views.views.RefsPage),
    url(r'^search', views.views.ResultsPage),
    url(r'^suggest', views.views.SuggestJson),
    url(r'^pathway$', views.views.DefinePathwayPage),
    url(r'^pathway/build_model', views.views.BuildPathwayModel),
    url(r'^pathway/results', views.views.PathwayResultPage),
    url(r'^robots\.txt$', TemplateView.as_view(template_name='robots.txt')),
    # Uncomment the next line to enable the admin:
    url(r'^admin/', include(admin.site.urls)),
]
