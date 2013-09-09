#!/usr/bin/python

import unittest
from util import inchi


class InchiTest(unittest.TestCase):
    
    def testAchiralInchi(self):
        # Test examples: input, expected output
        test_data = (('InChI=1/C3H7NO3/c4-2(1-5)3(6)7/h2,5H,1,4H2,(H,6,7)/t2-/m0/s1',
                      'InChI=1/C3H7NO3/c4-2(1-5)3(6)7/h2,5H,1,4H2,(H,6,7)'),
                     ('InChI=1/C3H7NO3/c4-2(1-5)3(6)7/h2,5H,1,4H2,(H,6,7)/t2-/m1/s1',
                      'InChI=1/C3H7NO3/c4-2(1-5)3(6)7/h2,5H,1,4H2,(H,6,7)'),
                     ('InChI=1/C6H12O6/c7-1-3-4(9)5(10)6(11,2-8)12-3/h3-5,7-11H,1-2H2/t3-,4-,5+,6-/m1/s1',
                      'InChI=1/C6H12O6/c7-1-3-4(9)5(10)6(11,2-8)12-3/h3-5,7-11H,1-2H2'),
                     ('InChI=1/C6H9N3O2/c7-5(6(10)11)1-4-2-8-3-9-4/h2-3,5H,1,7H2,(H,8,9)(H,10,11)/t5-/m1/s1/f/h8,10H',
                      'InChI=1/C6H9N3O2/c7-5(6(10)11)1-4-2-8-3-9-4/h2-3,5H,1,7H2,(H,8,9)(H,10,11)/f/h8,10H'),
                     ('InChI=1/C10H16O/c1-9(2)7-4-5-10(9,3)8(11)6-7/h7H,4-6H2,1-3H3/t7-,10+/m1/s1',
                      'InChI=1/C10H16O/c1-9(2)7-4-5-10(9,3)8(11)6-7/h7H,4-6H2,1-3H3'))
        
        for orig_inchi, achiral_inchi in test_data:
            self.assertEqual(achiral_inchi, inchi.AchiralInchi(orig_inchi))


if __name__ == '__main__':
    unittest.main()