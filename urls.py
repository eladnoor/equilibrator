import views.compound
import views.reaction
import views.views
import views.pathway
from django.conf.urls import include, url
from django.contrib import admin
from django.views.generic import TemplateView
admin.autodiscover()

urlpatterns = [
    url(r'^$', TemplateView.as_view(template_name='main.html')),
    url(r'^compound_image', views.compound.CompoundImage),
    url(r'^compound', views.compound.CompoundPage),
    url(r'^download', views.views.DownloadPage),
    url(r'^enzyme', views.views.EnzymePage),
    url(r'^reaction', views.reaction.ReactionPage),
    url(r'^graph_reaction', views.reaction.ReactionGraph),
    url(r'^data_refs', views.views.RefsPage),
    url(r'^search', views.views.ResultsPage),
    url(r'^suggest', views.views.SuggestJson),
    url(r'^pathway$', views.pathway.DefinePathwayPage),
    url(r'^pathway/build_model', views.pathway.BuildPathwayModel),
    url(r'^pathway/results', views.pathway.PathwayResultPage),
    url(r'^robots\.txt$', TemplateView.as_view(template_name='robots.txt')),
    # Uncomment the next line to enable the admin:
    url(r'^admin/', include(admin.site.urls)),
]
