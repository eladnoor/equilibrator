from django.shortcuts import render_to_response


def MainPage(request):
    """Renders the landing page."""
    return render_to_response('main.html', {})