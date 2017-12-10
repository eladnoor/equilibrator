#!/usr/bin/python3
import os
import logging
import itertools
import json
import re
from unittest import TestCase, main
from django.test import Client
import django

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "equilibrator.settings")
from equilibrator.settings import HAYSTACK_BACKEND

class GibbsTester(TestCase):

    def setUp(self):
        logging.getLogger().setLevel(logging.WARNING)
        django.setup()
        self.client = Client()

    def tearDown(self):
        pass

    def test_get_atom_bag(self):
        from gibbs import formula_parser
        from gibbs.models.compound import Compound

        # (formula, bag)
        test_data = (('C12H22O11', {'C': 12, 'H': 22, 'O': 11}),
                     ('C10H16N5O12P3S', {'C': 10, 'H': 16, 'N': 5,
                                         'O': 12, 'P': 3, 'S': 1}),
                     ('C10R11H16N5O12P3S', {'C': 10, 'H': 16, 'O': 12, 'N': 5, 'P': 3, 'S': 1, 'R': 11}),
                     ('C14H20O4(C5H8)n', {'H': 820, 'C': 514, 'O': 4}),
                     ('C34H32FeN4O4', {'H': 32, 'C': 34, 'Fe': 1, 'O': 4, 'N': 4}))
        
        parser = formula_parser.FormulaParser()
        compound = Compound(kegg_id='fake compound')
        for formula, expected_atom_bag in test_data:
            self.assertEqual(expected_atom_bag, parser.GetAtomBag(formula))

            compound.formula = formula
            if 'R' in expected_atom_bag:
                self.assertEqual(None, compound.GetAtomBag())
            else:
                self.assertEqual(expected_atom_bag, compound.GetAtomBag())        

    def test_transform(self):
        from gibbs.models.compound import Specie
        from gibbs.conditions import AqueousParams
        specie = Specie(number_of_hydrogens=2, number_of_mgs=0,
                        net_charge=-1, formation_energy=-10.0)
        
        # Test some hand-calculated numbers: (ph, ionic strength, result)
        test_data = (
             # Move pH, keep ionic strength constant.
             (6.0, 0.1, 59.071),
             (6.5, 0.1, 64.776),
             (7.0, 0.1, 70.481),
             (7.5, 0.1, 76.186),
             (8.0, 0.1, 81.891),
             # Move ionic strength, keep pH constant.
             (7.0, 0.0001, 69.898),
             (7.0, 0.001, 69.957),
             (7.0, 0.01, 70.121),
             (7.0, 0.05, 70.349),
             (7.0, 0.11, 70.501),
             (7.0, 0.15, 70.566),
             (7.0, 0.2, 70.629),
             # Move both.
             (6.0, 0.0001, 58.488),
             (6.5, 0.001, 64.252),
             (8.0, 0.0001, 81.308),
             (7.5, 0.15, 76.271),
             )

        for pH, ionic_strength, expected_transform in test_data:
            aq_params = AqueousParams(pH=pH, ionic_strength=ionic_strength)
            actual_transform = specie.Transform(aq_params)
            self.assertAlmostEqual(expected_transform, actual_transform, 3)

    def test_has_data(self):
        from gibbs.models.compound import Compound
        compound = Compound(kegg_id='fake compound')
        self.assertFalse(compound.HasData())
        
        compound.formula = 'C12H22O11'
        self.assertFalse(compound.HasData())

        compound.mass = 0.0
        self.assertFalse(compound.HasData())
        
        compound.mass = 14.5
        self.assertTrue(compound.HasData())
    
    def test_delta_g0(self):
        from gibbs.models.compound import Specie, SpeciesGroup, Compound
        from gibbs.conditions import AqueousParams
        from util import constants
        
        # Create a test compound.
        species = [Specie(number_of_hydrogens=12, net_charge=0,
                          formation_energy=-10.5),
                   Specie(number_of_hydrogens=11, net_charge=-1,
                          formation_energy=-12.1),
                   Specie(number_of_hydrogens=10, net_charge=-2,
                          formation_energy=-13.4)]
        species_group = SpeciesGroup()
        species_group._all_species = species

        compound = Compound(kegg_id='fake compound')
        compound._species_group_to_use = species_group
        
        # Format: ph, ionic strength, dG.
        test_data = ((6.5, 0.1, 361.094),
                     (7.0, 0.1, 389.619),
                     (7.5, 0.1, 418.143),
                     (7.0, 0.0001, 386.118),
                     (7.0, 0.001, 386.473),
                     (7.0, 0.2, 390.505))

        for pH, ionic_strength, expected_dg0 in test_data:
            aq_params = AqueousParams(pH=pH, ionic_strength=ionic_strength)
            actual_dg0 = compound.DeltaG0Prime(aq_params, phase=constants.DEFAULT_PHASE)
            self.assertAlmostEqual(expected_dg0, actual_dg0, 3,
                                   'ph: %f, i_s: %f, expected dG: %f, actual dG: %f' %
                                   (pH, ionic_strength, expected_dg0, actual_dg0))
    
    def test_hashable_reaction_string(self):
        """
            Ensure that hashable strings for different reactions are different
            and the same for equivalent reactions.
        """
        def reactant_list_to_hash(reactants):
            sparse = {r.compound.kegg_id : r.coeff for r in reactants}
            return StoredReaction.HashableReactionString(sparse)
        
        from gibbs.models.compound import Compound, Reactant
        from gibbs.models.reaction import StoredReaction
        
        compound_a = Compound(kegg_id='C00010')
        compound_b = Compound(kegg_id='C00009')
        compound_c = Compound(kegg_id='C00021')
        compound_d = Compound(kegg_id='C00032')
        compound_e = Compound(kegg_id='C00190')
        
        # 1 a + 1 b = 1 c + 1 d
        reactants_a = [Reactant(coeff=-1, compound=compound_a),
                       Reactant(coeff=-1, compound=compound_b),
                       Reactant(coeff=1, compound=compound_c),
                       Reactant(coeff=1, compound=compound_d)]
        sr_a = reactant_list_to_hash(reactants_a)
        
        # 1 a + 1 b = 1 c + 2 d
        reactants_b = [Reactant(coeff=-1, compound=compound_a),
                       Reactant(coeff=-1, compound=compound_b),
                       Reactant(coeff=1, compound=compound_c),
                       Reactant(coeff=2, compound=compound_d)]
        sr_b = reactant_list_to_hash(reactants_b)
        
        # 1 a + 1 b = 1 c + 1 e
        reactants_c = [Reactant(coeff=-1, compound=compound_a),
                       Reactant(coeff=-1, compound=compound_b),
                       Reactant(coeff=1, compound=compound_c),
                       Reactant(coeff=1, compound=compound_e)]
        sr_c = reactant_list_to_hash(reactants_c)
        
        # 3 a + 1 b = 2 c
        reactants_d = [Reactant(coeff=-3, compound=compound_a),
                       Reactant(coeff=-1, compound=compound_b),
                       Reactant(coeff=2, compound=compound_c)]
        sr_d = reactant_list_to_hash(reactants_d)

        # 1 a + 1 b = 1 d + 1 c
        reactants_a2 = [Reactant(coeff=-1, compound=compound_a),
                        Reactant(coeff=-1, compound=compound_b),
                        Reactant(coeff=1, compound=compound_d),
                        Reactant(coeff=1, compound=compound_c)]
        sr_a2 = reactant_list_to_hash(reactants_a2)
        
        # 1 d + 1 c = 1 a + 1 b 
        reactants_a3  = [Reactant(coeff=1, compound=compound_a),
                         Reactant(coeff=1, compound=compound_b),
                         Reactant(coeff=-1, compound=compound_d),
                         Reactant(coeff=-1, compound=compound_c)]
        sr_a3 = reactant_list_to_hash(reactants_a3)
        
        # 1 b + 3 a = 2 c
        reactants_d2 = [Reactant(coeff=-1, compound=compound_b),
                        Reactant(coeff=-3, compound=compound_a),
                        Reactant(coeff=2, compound=compound_c)]
        sr_d2 = reactant_list_to_hash(reactants_d2)
        
        # 2 c = 1 b + 3 a
        reactants_d3  = [Reactant(coeff=1, compound=compound_b),
                         Reactant(coeff=3, compound=compound_a),
                         Reactant(coeff=-2, compound=compound_c)]
        sr_d3 = reactant_list_to_hash(reactants_d3)

        # 1 b + 3 a = 2 c + 2 h+
        reactants_d4 = [Reactant(coeff=-1, compound=compound_b),
                        Reactant(coeff=-3, compound=compound_a),
                        Reactant(coeff=2, compound=compound_c)]
        sr_d4 = reactant_list_to_hash(reactants_d4)

        # make sure all the 4 reactions have different hashable strings
        for a, b in itertools.combinations((sr_a, sr_b, sr_c, sr_d), 2):
            self.assertNotEqual(a, b)
            
        # make sure all the 3 versions of reaction a have the same
        # hashable strings
        for a, b in itertools.combinations((sr_a, sr_a2, sr_a3), 2):
            self.assertEqual(a, b)

        # make sure all the 4 versions of reaction a have the same
        # hashable strings
        for a, b in itertools.combinations((sr_d, sr_d2, sr_d3, sr_d4), 2):
            self.assertEqual(a, b)

    def test_json_parsing(self):
        from gibbs.models.reaction import StoredReaction
        parsed_json = json.loads(
            """{ 
                "ECS": [ "4.1.1.39" ],
                "RID": "R00024",
                "names": [ "3-phospho-D-glycerate carboxy-lyase" ],          
                "reaction": [[-1, "C00001"], [-1, "C00011"],
                             [2.0, "C00197"], [-1, "C01182"]]  
                }
            """)
        expected_string = 'H2O + CO2 + d-ribulose 1,5-biphosphate = 2 3-Phospho-D-glycerate'
        
        rxn = StoredReaction.FromJson(parsed_json)
        self.assertEqual(expected_string, str(rxn))

    def testMinus(self):
        from gibbs.models.compound import Compound, CompoundWithCoeff
        from util import constants
        cmp = Compound(kegg_id='fake compound', formula='C10H16N5O12P3S')
        c_w_coeff = CompoundWithCoeff(coeff=4, compound=cmp, name='foo',
                                      phase=constants.DEFAULT_PHASE)
        minus_c = c_w_coeff.Minus()
        self.assertEqual(-4, minus_c.coeff)

    def test_glucose_using_kegg_id(self):
        response = self.client.post('/compound?compoundId=C00031&ph=8')
        self.assertEqual(response.status_code, 200)
        html = str(response.content).replace(r'\n', '')

        matches = re.findall('<td colspan="5">\s+<strong>([0-9\-\.]+)</strong>',
                             html)
            
        dgm, dg0 = map(float, matches)
        self.assertAlmostEqual(dgm, -378.4, 1)
        self.assertAlmostEqual(dg0, -361.3, 1)

    def test_atp_search_results(self):
        response = self.client.post('/search?query=ATP')
        self.assertEqual(response.status_code, 200)
        
        html = str(response.content).replace(r'\n', '')

        # match the compound result
        matches = re.findall('<tr class="infoTableHeader">\s+<th colspan="100%">(\w+)</th>',
                             html)
        self.assertIn('ATP', matches)
        self.assertIn('dATP', matches)
        self.assertIn('ATPA', matches)
        
        if HAYSTACK_BACKEND == 'solr':
            # the search doesn't work very well with "simple" and these
            # results are not in the matched list.
            self.assertIn('BzATP', matches)

            # match the enzyme result
            matches = re.findall('<tr class="infoTableHeader">\s+<th colspan="100%">\s+<a href="[^"]+">([\w\-]+)</a>',
                                 html)
            self.assertIn('AtPAO3', matches)
            self.assertIn('AtPAO3', matches)
            self.assertIn('ATPase', matches)


if __name__ == "__main__":
    main()
