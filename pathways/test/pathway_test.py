from django.test import TestCase

import pathways

class PathwayTest(TestCase):

    def test_from_file(self):
        path = pathways.Pathway.from_filename('pathways/test/simple_pathway.csv')
        path.print_reactions()