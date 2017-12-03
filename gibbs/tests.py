from django.test import TestCase

import unittest
import json
import itertools
from gibbs import models_test
from gibbs import reaction_test
from gibbs import formula_parser
from gibbs.reaction import Reaction
from gibbs import models
from gibbs import reaction

# Create your tests here.
def gibbs_test():
    test_modules = (models_test,
                    reaction_test)
    
    modules_str = ', '.join(m.__name__ for m in test_modules)
    print('Running test suites from modules %s' % modules_str)
    
    suites = [m.Suite() for m in test_modules]
    alltests = unittest.TestSuite(suites)
    
    runner = unittest.TextTestRunner()
    runner.run(alltests)
    
class FormulaParserTest(unittest.TestCase):
    
    def testGetAtomBag(self):
        # (formula, bag)
        test_data = ((None, None),
                     ('C12H22O11', {'C': 12, 'H': 22, 'O': 11}),
                     ('C10H16N5O12P3S', {'C': 10, 'H': 16, 'N': 5,
                                         'O': 12, 'P': 3, 'S': 1}),
                     ('C10R11H16N5O12P3S', {'C': 10, 'H': 16, 'O': 12, 'N': 5, 'P': 3, 'S': 1, 'R': 11}),
                     ('C14H20O4(C5H8)n', {'H': 820, 'C': 514, 'O': 4}),
                     ('C34H32FeN4O4', {'H': 32, 'C': 34, 'Fe': 1, 'O': 4, 'N': 4}))
        
        parser = formula_parser.FormulaParser()
        for formula, expected_bag in test_data:
            #parser.GetAtomBag(formula)
            self.assertEqual(expected_bag, parser.GetAtomBag(formula))
  
 
class SpeciesFormationEnergyTest(unittest.TestCase):
    """Tests for SpeciesFormationEnergy"""
        
    def testTransform(self):
        specie = models.Specie(number_of_hydrogens=2, number_of_mgs=0,
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

        for ph, ionic_strength, expected_transform in test_data:
            actual_transform = specie.Transform(pH=ph, ionic_strength=ionic_strength)
            self.assertAlmostEqual(expected_transform, actual_transform, 3)


class CompoundTest(unittest.TestCase):
    
    def testHasData(self):
        compound = models.Compound(kegg_id='fake compound')
        self.assertFalse(compound.HasData())
        
        compound.formula = 'C12H22O11'
        self.assertFalse(compound.HasData())

        compound.mass = 0.0
        self.assertFalse(compound.HasData())
        
        compound.mass = 14.5
        self.assertTrue(compound.HasData())
    
    def testGetAtomBag(self):
        compound = models.Compound(kegg_id='fake compound')
        self.assertEqual(None, compound.GetAtomBag())
        
        compound.formula = 'C12H22O11'
        expected_atom_bag = {'C': 12, 'H': 22, 'O': 11}
        self.assertEqual(expected_atom_bag, compound.GetAtomBag())
        
        compound.formula = 'C10H16N5O12P3S'
        expected_atom_bag = {'C': 10, 'H': 16, 'N': 5,
                             'O': 12, 'P': 3, 'S': 1}
        self.assertEqual(expected_atom_bag, compound.GetAtomBag())

        # Contains an R group.
        compound.formula = 'C10R11H16N5O12P3S'
        self.assertEqual(None, compound.GetAtomBag())

    def testDeltaG0(self):
        # Create a test compound.
        species = [models.Specie(number_of_hydrogens=12, net_charge=0,
                                 formation_energy=-10.5),
                   models.Specie(number_of_hydrogens=11, net_charge=-1,
                                 formation_energy=-12.1),
                   models.Specie(number_of_hydrogens=10, net_charge=-2,
                                 formation_energy=-13.4)]
        species_group = models.SpeciesGroup()
        species_group._all_species = species

        compound = models.Compound(kegg_id='fake compound')
        compound._species_group_to_use = species_group
        
        # Format: ph, ionic strength, dG.
        test_data = ((6.5, 0.1, 361.094),
                     (7.0, 0.1, 389.619),
                     (7.5, 0.1, 418.143),
                     (7.0, 0.0001, 386.118),
                     (7.0, 0.001, 386.473),
                     (7.0, 0.2, 390.505))

        for ph, i_s, expected_dg0 in test_data:
            actual_dg0 = compound.DeltaG0Prime(pH=ph, ionic_strength=i_s)
            self.assertAlmostEqual(expected_dg0, actual_dg0, 3,
                                   'ph: %f, i_s: %f, expected dG: %f, actual dG: %f' %
                                   (ph, i_s, expected_dg0, actual_dg0))

class StoredReactionTest(unittest.TestCase):
    """Tests for StoredReaction"""
    
    def testHashableReactionString(self):
        """
            Ensure that hashable strings for different reactions are different
            and the same for equivalent reactions.
        """
        
        compound_a = models.Compound(kegg_id='C00010')
        compound_b = models.Compound(kegg_id='C00009')
        compound_c = models.Compound(kegg_id='C00021')
        compound_d = models.Compound(kegg_id='C00032')
        compound_e = models.Compound(kegg_id='C00190')
        
        # 1 a + 1 b = 1 c + 1 d
        reactants_a = [models.Reactant(coeff=-1, compound=compound_a),
                       models.Reactant(coeff=-1, compound=compound_b),
                       models.Reactant(coeff=1, compound=compound_c),
                       models.Reactant(coeff=1, compound=compound_d)]
        sr_a = models.StoredReaction.HashableReactionString(reactants_a)
        
        # 1 a + 1 b = 1 c + 2 d
        reactants_b = [models.Reactant(coeff=-1, compound=compound_a),
                       models.Reactant(coeff=-1, compound=compound_b),
                       models.Reactant(coeff=1, compound=compound_c),
                       models.Reactant(coeff=2, compound=compound_d)]
        sr_b = models.StoredReaction.HashableReactionString(reactants_b)
        
        # 1 a + 1 b = 1 c + 1 e
        reactants_c = [models.Reactant(coeff=-1, compound=compound_a),
                       models.Reactant(coeff=-1, compound=compound_b),
                       models.Reactant(coeff=1, compound=compound_c),
                       models.Reactant(coeff=1, compound=compound_e)]
        sr_c = models.StoredReaction.HashableReactionString(reactants_c)
        
        # 3 a + 1 b = 2 c
        reactants_d = [models.Reactant(coeff=-3, compound=compound_a),
                       models.Reactant(coeff=-1, compound=compound_b),
                       models.Reactant(coeff=2, compound=compound_c)]
        sr_d = models.StoredReaction.HashableReactionString(reactants_d)

        # 1 a + 1 b = 1 d + 1 c
        reactants_a2 = [models.Reactant(coeff=-1, compound=compound_a),
                        models.Reactant(coeff=-1, compound=compound_b),
                        models.Reactant(coeff=1, compound=compound_d),
                        models.Reactant(coeff=1, compound=compound_c)]
        sr_a2 = models.StoredReaction.HashableReactionString(reactants_a2)
        
        # 1 d + 1 c = 1 a + 1 b 
        reactants_a3  = [models.Reactant(coeff=1, compound=compound_a),
                         models.Reactant(coeff=1, compound=compound_b),
                         models.Reactant(coeff=-1, compound=compound_d),
                         models.Reactant(coeff=-1, compound=compound_c)]
        sr_a3 = models.StoredReaction.HashableReactionString(reactants_a3)
        
        # 1 b + 3 a = 2 c
        reactants_d2 = [models.Reactant(coeff=-1, compound=compound_b),
                        models.Reactant(coeff=-3, compound=compound_a),
                        models.Reactant(coeff=2, compound=compound_c)]
        sr_d2 = models.StoredReaction.HashableReactionString(reactants_d2)
        
        # 2 c = 1 b + 3 a
        reactants_d3  = [models.Reactant(coeff=1, compound=compound_b),
                         models.Reactant(coeff=3, compound=compound_a),
                         models.Reactant(coeff=-2, compound=compound_c)]
        sr_d3 = models.StoredReaction.HashableReactionString(reactants_d3)

        # 1 b + 3 a = 2 c + 2 h+
        reactants_d4 = [models.Reactant(coeff=-1, compound=compound_b),
                        models.Reactant(coeff=-3, compound=compound_a),
                        models.Reactant(coeff=2, compound=compound_c)]
        sr_d4 = models.StoredReaction.HashableReactionString(reactants_d4)

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
        

class ReactionTest(unittest.TestCase):

    def testJsonParsing(self):
        parsed_json = json.loads(
            """{ 
                "ECS": [ "4.1.1.39" ],
                "RID": "R00024",
                "names": [ "3-phospho-D-glycerate carboxy-lyase" ],          
                "reaction": [[-1, "C00001"], [-1, "C00011"],
                             [2.0, "C00197"], [-1, "C01182"]]  
                }
            """)
        rxn = Reaction.FromJson(parsed_json)
        print(str(rxn))

class CompoundWithCoeffTest(unittest.TestCase):
    
    def setUp(self):
        self.compound = models.Compound(kegg_id='fake compound',
                                        formula='C10H16N5O12P3S')
    
    def testMinus(self):
        c_w_coeff = reaction.CompoundWithCoeff(coeff=4, compound=self.compound,
                                               name='foo')
        minus_c = c_w_coeff.Minus()
        self.assertEquals(-4, minus_c.coeff)
        
    def testMicromolarConcentration(self):
        c_w_coeff = reaction.CompoundWithCoeff(coeff=4, compound=self.compound,
                                               name='foo', concentration=0.5)
        self.assertAlmostEqual(5e5, c_w_coeff.micromolar_concentration, 3)

def Suite():
    suites = (unittest.makeSuite(SpeciesFormationEnergyTest, 'test'),
              unittest.makeSuite(CompoundTest, 'test'),
              unittest.makeSuite(StoredReactionTest, 'test'),
              unittest.makeSuite(ReactionTest, 'test'),
              unittest.makeSuite(FormulaParserTest, 'test'),
              unittest.makeSuite(CompoundWithCoeffTest, 'test')
              )
    return unittest.TestSuite(suites)
