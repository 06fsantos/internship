'''
Created on 9 Apr 2018

@author: filipe
'''
import unittest
import codon_frequency_table


class Test(unittest.TestCase):

    
    def test_update_dict(self):
        file = 'individual_symbol_for_test.xlsx'
        sheet = 'Sheet1'
        codon_dict = {
                    'AUA':0, 'AUC':0, 'AUU':0, 'AUG':0,
                    'ACA':0, 'ACC':0, 'ACG':0, 'ACU':0, 
                    'AAC':0, 'AAU':0, 'AAA':0, 'AAG':0,
                    'AGC':0, 'AGU':0, 'AGA':0, 'AGG':0,
                    'CUA':0, 'CUC':0, 'CUG':0, 'CUU':0,
                    'CCA':0, 'CCC':0, 'CCG':0, 'CCU':0,
                    'CAC':0, 'CAU':0, 'CAA':0, 'CAG':0,
                    'CGA':0, 'CGC':0, 'CGG':0, 'CGU':0,
                    'GUA':0, 'GUC':0, 'GUG':0, 'GUU':0,
                    'GCA':0, 'GCC':0, 'GCG':0, 'GCU':0,
                    'GAC':0, 'GAU':0, 'GAA':0, 'GAG':0,
                    'GGA':0, 'GGC':0, 'GGG':0, 'GGU':0,
                    'UCA':0, 'UCC':0, 'UCG':0, 'UCU':0,
                    'UUC':0, 'UUU':0, 'UUA':0, 'UUG':0,
                    'UAC':0, 'UAU':0, 'UGC':0, 'UGU':0, 
                    'UGG':0}
        
        new_table = {'AUA': 0.0, 'AUC': 0.001, 'AUU': 0.0, 'AUG': 0.001,
                    'ACA': 0.001, 'ACC': 0.0, 'ACG': 0.0, 'ACU': 0.002, 
                    'AAC': 0.002, 'AAU': 0.0, 'AAA': 0.0, 'AAG': 0.007,
                    'AGC': 0.004, 'AGU': 0.001, 'AGA': 0.0, 'AGG': 0.0,
                    'CUA': 0.0, 'CUC': 0.002, 'CUG': 0.005, 'CUU': 0.0, 
                    'CCA': 0.002, 'CCC': 0.003, 'CCG': 0.001, 'CCU': 0.0, 
                    'CAC': 0.003, 'CAU': 0.001, 'CAA': 0.003, 'CAG': 0.003, 
                    'CGA': 0.0, 'CGC': 0.001, 'CGG': 0.0, 'CGU': 0.001, 
                    'GUA': 0.0, 'GUC': 0.0, 'GUG': 0.003, 'GUU': 0.0, 
                    'GCA': 0.0, 'GCC': 0.005, 'GCG': 0.0, 'GCU': 0.0, 
                    'GAC': 0.002, 'GAU': 0.001, 'GAA': 0.008, 'GAG': 0.019, 
                    'GGA': 0.0, 'GGC': 0.005, 'GGG': 0.001, 'GGU': 0.002, 
                    'UCA': 0.002, 'UCC': 0.005, 'UCG': 0.0, 'UCU': 0.0, 
                    'UUC': 0.0, 'UUU': 0.001, 'UUA': 0.0, 'UUG': 0.0, 
                    'UAC': 0.0, 'UAU': 0.0, 'UGC': 0.0, 'UGU': 0.002, 
                    'UGG': 0.0}
        
        self.assertEqual(codon_frequency_table.update_dict(file, sheet, codon_dict), new_table, msg = 'Error: incorrect frequency count')


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_update_dict']
    unittest.main()