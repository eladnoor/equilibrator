#!/usr/bin/python

import unittest
from util import kegg


class KeggTest(unittest.TestCase):
    
    def testKeggIdFromInt(self):
        # Test examples: input, expected output
        test_data = ((1, 'C00001'), (80, 'C00080'), (1985, 'C01985'))
        
        for id, str_id in test_data:
            self.assertEqual(str_id, kegg.KeggIdFromInt(id))
        

if __name__ == '__main__':
    unittest.main()