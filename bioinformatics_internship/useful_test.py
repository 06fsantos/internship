'''
Created on 28 Feb 2018

@author: filipe
'''
import unittest
import useful

class Test(unittest.TestCase):


    def test_pull_fasta_sequence(self):
        gene_id = 'NM_029459'
        self.assertEqual(useful.pull_fasta_sequence(gene_id), 'AGGGACCGTTGAGGGGCAGCTTCCACCAAAGACTATGGCACGCCCACCACCTCGAACTCCTCTCCAGAAA\
TGAACGACTAACACTGCTGAGGGAGTGGAAAAGTTAAAAAAAATAAGAAGAAAGGAAAAAAAAGAAAGAA\
AAGAAAAAAAACGAAAAGAAAATCTCCAAGTCCGCCCACTTCTTCATTACATCCTTTCACTCCCTCTTCC\
AAGATTTATATTCCTGAACCATTCGTGGTGGCAATTCCTGGGCTGTGAGAGGAATTTCGAGGTCTGCGTC\
AACTGGGATAACGAAAGTGGACATTATTTCAAATATTTCATTGAATTTGATCAGCGTTTTTCCATAGTCT\
CATCCAGAGAGATAGATCTTCACTGGATTCACAACTCAGACACATTTGAAGATTCTTGTAGAGCATCTGT\
GAGAGGAAGGAGGCTGCTGCAACCTAAGGCCTTTGTGGGTCTGGAACTCAGGAATCTCAGTTTCTGCAAT\
CTTGACCCTTACTGAAGTGCACCGGTTCCAGGAACCCATCAAATTTGAGTGATATTTGAAAGCCCTTTGC\
TACGGTGAGATTAATGAAGAGCTGTCCTGTATTCTTCAAGAAGCTGGTAATATTTATTTCAGTCAGCAAG\
TTATTCTAAACAAGAACAGTGTCTGAGTGGCAAGTTATTCTATCCAAGGACGGTGTCTGAGTGTGTACCG\
GCTAATAGTAAAGTTCCCTAAACTAGGTTTATGATGATGGATGAAAGAGACCCATCCTCGCTTTTGGATC\
TGGCTATACAGAGTCTACTAAGTAATGAGCTTGTAGCAATTCATTCTCTGGGGGAGATCCCAAGGGAGCT\
TTTTGTTCCATTGTTCTCTGCTGCCTTCACGGGAGGATATAGGAAGATACTGACTTCAATGGTGAAGATT\
TGGCCTTTTACCTGTCTCCACATTGGAACATTAAGTGTACAGGAACCCCAGCGTGAACTCCTGAAAGCCA\
TGGTTGAGAGTCTTCAGTTTCTTCCTGCCCAGGACTGTTCTTCTGGGGGCCCTAAGTTGAGGATCCTAGA\
TGTAAGGCAGGGTGTTGACTGCAAGACAACATGCCCTGATTTTGGTGCCAGATCTCCAACTTGTTTTCAT\
GGTTGTACTCACTCTGTACACTCTATTCTGAAGTTAGAAAGCCAGTACAGCATTGTAGATCTAAAGCCCG\
AGAGTCAGTCTGCAATCCAGCCTATGGAACTACTAGTAGACCTTTCCCTTGATGGTACCTTGAGAGAAAG\
GGAATTTTTTGCTTTGCTTCTGAATAAAGTACAGCAGAGCTCAGGGTCTTTGCACCTCTGCTGCCGAGAT\
CTACAAATTGATAGATTTTCTTATGCCAAAAACGCTCTGAAGTTCCTCGATCTAACTTGCATTCAGAACC\
TGACAGTTGATCAGGCTTCACTGAGTGAAGTCACCACTCTTCTGGCTCGCATGATCTATCTGGACAGCCT\
GAGTCTCTCTAAAATCACTTATAGATCTTTGCATGGGAAAGTCTTCCGAGTGTTCCTCAACTATCTTGGG\
CGGATGAACTGCCTGAAAGAGCTCAACCTGTCTTCCTTTAGCCTCACAGACCATCTGGATAGCCTCCTCA\
GAGCCTTACCACCTAATTTGGATTTCTTGTATCTGCCGTTCTGTGAAATTTCTTACAGAGATCTCAAATT\
TCTATCCCAGAGTGCTCAGGCCACCCACTTAAAGCTGTTGAATCTTAGTAACAACCCAATGTATTGGGAT\
GATTGTGGGCCTTTTCAGACTCTTTTGCAGAAGCTCTCAGATACCTTGCAGCATCTGGCCATAAACCATT\
GCCATTTAACAGATGCCATACTCTCTGCTATTCTGCCAGCACTATCTAAGTGTTCCCATCTCCGTGTGAT\
TAGCTTTGTCTCTAACCCCATTTCAATGCCTATGCTCCTGAAAATTCTTCATTACTTAACACCTTTGATG\
GAGCTGAAATACGTGATTTACCCTATCCCTATACATTGCTATGAACAATGGCAATTTCATGGCAGATTAG\
ACCGGCAGAAGCTCACCGATGTCCAAGCACAACTGAAGGCAATGCTACAAGCAGCAAAAAGGAGTGACAT\
GAACTGGATCACTTATTCTCAGTAAACTTCCAAGTTTAACTCCATCTCAAGCTCCAAATTTGACCTGTTA\
TCTGTTCAATGTTCTTTTCTCGAGCTTCAAGAATCTGATGTAAGAATTCGTACGTTATAGACGATTAAAG\
TTAGAAACTGATCAAAAACATTAACTC', msg = 'Error: output not equal to sequence')
        self.assertNotEqual(useful.pull_fasta_sequence(gene_id), 'AGGGACCGTTGAGGGGCAGCTTCCACCAAAGACTATGGCACGCCCACCACCTCGAACTCCTCTCCAGAAA\
TGAACGACTAACACTGCTGAGGGAGTGGAAAAGTTAAAAAAAATAAGAAGAAAGGAAAAAAAAGAAAGAA\
AAGAAAAAAAACGAAAAGAAAATCTCCAAGTCCGCCCACTTCTTCATTACATCCTTTCACTCCCTCTTCC\
AAGATTTATATTCCTGAACCATTCGTGGTGGCAATTCCTGGGCTGTGAGAGGAATTTCGAGGTCTGCGTC\
AACTGGGATAACGAAAGTGGACATTATTTCAAATATTTCATTGAATTTGATCAGCGTTTTTCCATAGTCT\
CATCCAGAGAGATAGATCTTCACTGGATTCACAACTCAGACACATTTGAAGATTCTTGTAGAGCATCTGT\
GAGAGGAAGGAGGCTGCTGCAACCTAAGGCCTTTGTGGGTCTGGAACTCAGGAATCTCAGTTTCTGCAAT\
CTTGACCCTTACTGAAGTGCACCGGTTCCAGGAACCCATCAAATTTGAGTGATATTTGAAAGCCCTTTGC\
TACGGTGAGATTAATGAAGAGCTGTCCTGTATTCTTCAAGAAGCTGGTAATATTTATTTCAGTCAGCAAG\
TTATTCTAAACAAGAACAGTGTCTGAGTGGCAAGTTATTCTATCCAAGGACGGTGTCTGAGTGTGTACCG\
GCTAATAGTAAAGTTCCCTAAACTAGGTTTATGATGATGGATGAAAGAGACCCATCCTCGCTTTTGGATC\
TGGCTATACAGAGTCTACTAAGTAATGAGCTTGTAGCAATTCATTCTCTGGGGGAGATCCCAAGGGAGCT\
TTTTGTTCCATTGTTCTCTGCTGCCTTCACGGGAGGATATAGGAAGATACTGACTTCAATGGTGAAGATT\
TGGCCTTTTACCTGTCTCCACATTGGAACATTAAGTGTACAGGAACCCCAGCGTGAACTCCTGAAAGCCA\
TGGTTGAGAGTCTTCAGTTTCTTCCTGCCCAGGACTGTTCTTCTGGGGGCCCTAAGTTGAGGATCCTAGA\
TGTAAGGCAGGGTGTTGACTGCAAGACAACATGCCCTGATTTTGGTGCCAGATCTCCAACTTGTTTTCAT\
GGTTGTACTCACTCTGTACACTCTATTCTGAAGTTAGAAAGCCAGTACAGCATTGTAGATCTAAAGCCCG\
AGAGTCAGTCTGCAATCCAGCCTATGGAACTACTAGTAGACCTTTCCCTTGATGGTACCTTGAGAGAAAG\
GGAATTTTTTGCTTTGCTTCTGAATAAAGTACAGCAGAGCTCAGGGTCTTTGCACCTCTGCTGCCGAGAT\
CTACAAATTGATAGATTTTCTTATGCCAAAAACGCTCTGAAGTTCCTCGATCTAACTTGCATTCAGAACC\
TGACAGTTGATCAGGCTTCACTGAGTGAAGTCACCACTCTTCTGGCTCGCATGATCTATCTGGACAGCCT\
GAGTCTCTCTAAAATCACTTATAGATCTTTGCATGGGAAAGTCTTCCGAGTGTTCCTCAACTATCTTGGG\
CGGATGAACTGCCTGAAAGAGCTCAACCTGTCTTCCTTTAGCCTCACAGACCATCTGGATAGCCTCCTCA\
GAGCCTTACCACCTAATTTGGATTTCTTGTATCTGCCGTTCTGTGAAATTTCTTACAGAGATCTCAAATT\
TCTATCCCAGAGTGCTCAGGCCACCCACTTAAAGCTGTTGAATCTTAGTAACAACCCAATGTATTGGGAT\
GATTGTGGGCCTTTTCAGACTCTTTTGCAGAAGCTCTCAGATACCTTGCAGCATCTGGCCATAAACCATT\
GCCATTTAACAGATGCCATACTCTCTGCTATTCTGCCAGCACTATCTAAGTGTTCCCATCTCCGTGTGAT\
TAGCTTTGTCTCTAACCCCATTTCAATGCCTATGCTCCTGAAAATTCTTCATTACTTAACACCTTTGATG\
GAGCTGAAATACGTGATTTACCCTATCCCTATACATTGCTATGAACAATGGCAATTTCATGGCAGATTAG\
ACCGGCAGAAGCTCACCGATGTCCAAGCACAACTGAAGGCAATGCTACAAGCAGCAAAAAGGAGTGACAT\
GAACTGGATCACTTATTCTCAGTAAACTTCCAAGTTTAACTCCATCTCAAGCTCCAAATTTGACCTGTTA\
TCTGTTCAATGTTCTTTTCTCGAGCTTCAAGAATCTGATGTAAGAATTCGTACGTTATAGACGATTAAAG\
TTAGAAACTGATCAAAAACATTAAAAA', msg = 'Error: does not recognise different sequences')
    
    def test_count_codons(self):
        sequence = 'AAATGACGACGATTACG'
        codon = 'ACG'
        self.assertEqual(useful.count_codons(sequence, codon), 3, msg = 'Error: incorrect count')
        self.assertNotEqual(useful.count_codons(sequence, codon), 4, msg = 'Error: counting too many codons' )
        ######### check if only counts nucleotides in frame ########
        sequence = 'ATCATGAACGGCAACATAACG' 
        codon = 'AAC'
        self.assertEqual(useful.count_codons(sequence, codon), 2, 'Error: counting nucleotides out of frame')  
        ######## check if stop codon breaks the loop 
        sequence = 'ATGAACAACTAAAAC'
        self.assertEqual(useful.count_codons(sequence, codon), 2, 'Error: Stop codon not recognised')

    def test_clean_seq(self):
        sequence = 'aacggttaa'
        self.assertEqual(useful.clean_seq(sequence), 'AACGGTTAA', msg = 'Error: does not utilise .upper()')
        sequence = 'aaggttddaatt'
        self.assertEqual(useful.clean_seq(sequence), 'AAGGTTAATT', msg = 'Error: does not remove non nucleotide letters')
        self.assertEqual(useful.clean_seq(' '), '', msg = 'Error: does not recognises spaces to skip')

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test']
    unittest.main()