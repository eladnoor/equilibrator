import logging

from django.http import Http404
from django.shortcuts import render_to_response
from gibbs.forms import EnzymeForm
from gibbs import models
from util import django_utils


def EnzymePage(request):    
    """Renders a page for a particular enzyme."""
    form = EnzymeForm(request.GET)
    if not form.is_valid():
        logging.error(form.errors)
        raise Http404

    enz = models.Enzyme.objects.get(ec=form.cleaned_ec)
    template_data = {'is_superuser': django_utils.IsSuperUser(request),
                     'enzyme': enz}
    return render_to_response('enzyme_page.html', template_data)