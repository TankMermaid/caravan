#!/usr/bin/env python

'''
unit tests for ???
'''

import unittest
from caravan.test import fake_fh
from caravan import barcodes

class TestParseBarcode(unittest.TestCase):
    def test_correct(self):
        '''make sure the function provides the correct output'''
        self.assertEqual(barcodes.MappedRecords.parse_barcode('@any_set_of_chars#ACGT/1'), 'ACGT')
    
    def test_raise(self):
        '''make sure the function raises an error when appropriate'''
        pass


if __name__ == '__main__':
    unittest.main(verbosity=2)
