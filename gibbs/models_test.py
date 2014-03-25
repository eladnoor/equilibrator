#!/usr/bin/python

import itertools
import unittest
from util import django_utils

# NOTE(flamholz): This is crappy. We're using the real database for
# a unit test. I wish I knew of a better way.
django_utils.SetupDjango()

from gibbs import models

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
            actual_dg0 = compound.DeltaG0(pH=ph, ionic_strength=i_s)
            self.assertAlmostEqual(expected_dg0, actual_dg0, 3,
                                   'ph: %f, i_s: %f, expected dG: %f, actual dG: %f' %
                                   (ph, i_s, expected_dg0, actual_dg0))

class StoredReactionTest(unittest.TestCase):
    """Tests for StoredReaction"""
    
    def testHashableReactionString(self):
        """Ensure that hashable strings for different reactions are different."""
        compound_a = models.Compound(kegg_id='C00010')
        compound_b = models.Compound(kegg_id='C00009')
        compound_c = models.Compound(kegg_id='C00021')
        compound_d = models.Compound(kegg_id='C00032')
        compound_e = models.Compound(kegg_id='C00190')
        hydrogen = models.Compound(kegg_id='C00080')
        
        # 1 a + 1 b = 1 c + 1 d
        reactants_a = [models.Reactant(coeff=1, compound=compound_a),
                       models.Reactant(coeff=1, compound=compound_b)]
        products_a  = [models.Reactant(coeff=1, compound=compound_c),
                       models.Reactant(coeff=1, compound=compound_d)]
        hashable_a  = models.StoredReaction.HashableReactionString(reactants_a,
                                                                   products_a)
        
        
        # 1 a + 1 b = 1 c + 2 d
        reactants_b = [models.Reactant(coeff=1, compound=compound_a),
                       models.Reactant(coeff=1, compound=compound_b)]
        products_b  = [models.Reactant(coeff=1, compound=compound_c),
                       models.Reactant(coeff=2, compound=compound_d)]
        hashable_b  = models.StoredReaction.HashableReactionString(reactants_b,
                                                                   products_b)
        
        # 1 a + 1 b = 1 c + 1 e
        reactants_c = [models.Reactant(coeff=1, compound=compound_a),
                       models.Reactant(coeff=1, compound=compound_b)]
        products_c  = [models.Reactant(coeff=1, compound=compound_c),
                       models.Reactant(coeff=1, compound=compound_e)]
        hashable_c  = models.StoredReaction.HashableReactionString(reactants_c,
                                                                   products_c)
        
        # 3 a + 1 b = 2 c
        reactants_d = [models.Reactant(coeff=3, compound=compound_a),
                       models.Reactant(coeff=1, compound=compound_b)]
        products_d  = [models.Reactant(coeff=2, compound=compound_c)]
        hashable_d  = models.StoredReaction.HashableReactionString(reactants_d,
                                                                   products_d)
        
        all_hashables = (hashable_a, hashable_b, hashable_c, hashable_d)
        for a, b in itertools.combinations(all_hashables, 2):
            self.assertNotEqual(a, b)
            
        # 1 a + 1 b = 1 d + 1 c
        reactants_a2 = [models.Reactant(coeff=1, compound=compound_a),
                        models.Reactant(coeff=1, compound=compound_b)]
        products_a2  = [models.Reactant(coeff=1, compound=compound_d),
                        models.Reactant(coeff=1, compound=compound_c)]
        hashable_a2  = models.StoredReaction.HashableReactionString(reactants_a2,
                                                                    products_a2)
        self.assertEqual(hashable_a2, hashable_a)
        
        # 1 d + 1 c = 1 a + 1 b 
        products_a3  = [models.Reactant(coeff=1, compound=compound_a),
                        models.Reactant(coeff=1, compound=compound_b)]
        reactants_a3 = [models.Reactant(coeff=1, compound=compound_d),
                        models.Reactant(coeff=1, compound=compound_c)]
        hashable_a3  = models.StoredReaction.HashableReactionString(reactants_a3,
                                                                    products_a3)
        self.assertEqual(hashable_a3, hashable_a)
        
        # 1 b + 3 a = 2 c
        reactants_d2 = [models.Reactant(coeff=1, compound=compound_b),
                        models.Reactant(coeff=3, compound=compound_a)]
        products_d2  = [models.Reactant(coeff=2, compound=compound_c)]
        hashable_d2  = models.StoredReaction.HashableReactionString(reactants_d2,
                                                                    products_d2)
        self.assertEqual(hashable_d2, hashable_d)
        
        # 2 c = 1 b + 3 a
        products_d3  = [models.Reactant(coeff=1, compound=compound_b),
                        models.Reactant(coeff=3, compound=compound_a)]
        reactants_d3 = [models.Reactant(coeff=2, compound=compound_c)]
        hashable_d3  = models.StoredReaction.HashableReactionString(reactants_d3,
                                                                    products_d3)
        self.assertEqual(hashable_d3, hashable_d)

        # 1 b + 3 a = 2 c + 2 h+
        reactants_d4 = [models.Reactant(coeff=1, compound=compound_b),
                        models.Reactant(coeff=3, compound=compound_a)]
        products_d4  = [models.Reactant(coeff=2, compound=compound_c),
                        models.Reactant(coeff=2, compound=hydrogen)]
        hashable_d4  = models.StoredReaction.HashableReactionString(reactants_d4,
                                                                    products_d4)
        self.assertEqual(hashable_d4, hashable_d)
        

def Suite():
    suites = (unittest.makeSuite(SpeciesFormationEnergyTest, 'test'),
              unittest.makeSuite(CompoundTest, 'test'),
              unittest.makeSuite(StoredReactionTest, 'test'))
    return unittest.TestSuite(suites)
    

if __name__ == '__main__':
    unittest.main()