'''
Created on 18 May 2018

@author: filipe
'''
import useful 
import pandas as pd 
from Bio import Seq

def count_gc(file, sheet, index_column):
    
    gc_count = 0
    total_count = 0 
    
    df = pd.read_excel(file, sheetname = sheet , index_col = index_column)
    symbol = df['Symbol'].copy()
    
    for chromosome_id in symbol:
        seq = useful.pull_fasta_sequence(chromosome_id)
        seq = useful.clean_seq(seq)
        seq = Seq.Seq(seq, Seq.Alphabet.generic_dna)
        
        for nucleotide in seq:
            if nucleotide == 'G' or nucleotide == 'C':
                gc_count += 1
                total_count += 1
            else:
                total_count += 1
                
    gc_percentage = (gc_count / total_count) * 100 
    return gc_percentage  
    
if __name__ == '__main__':
    print ('Cell Lines - Mock vs Ala up: {}'.format(count_gc('cell_lines_clean.xlsx', 'Mock vs Ala up', None)))
    print ('Cell Lines - Mock vs Ala down: {}'.format(count_gc('cell_lines_clean.xlsx', 'Mock vs Ala down', None)))
    print ('Cell Lines - Mock vs Wt up: {}'.format(count_gc('cell_lines_clean.xlsx', 'Mock vs Wt up', None)))
    print ('Cell Lines - Mock vs Wt down: {}'.format(count_gc('cell_lines_clean.xlsx', 'Mock vs Wt down', None)))