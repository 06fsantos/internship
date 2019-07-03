'''
Created on 26 Jun 2018

@author: filipe
'''

import useful
import pandas as pd 
from Bio import Seq
from collections import Counter

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


def get_rscu_value(codon_count_dict):

    synonymous_codon_dict = { 
    'CYS': ['UGU', 'UGC'], 
    'ASP': ['GAU', 'GAC'], 
    'SER': ['UCU', 'UCG', 'UCA', 'UCC', 'AGC', 'AGU'], 
    'GLN': ['CAA', 'CAG'], 
    'MET': ['AUG'], 
    'ASN': ['AAC', 'AAU'], 
    'PRO': ['CCU', 'CCG', 'CCA', 'CCC'], 
    'LYS': ['AAG', 'AAA'],  
    'THR': ['ACC', 'ACA', 'ACG', 'ACU'], 
    'PHE': ['UUU', 'UUC'], 
    'ALA': ['GCA', 'GCC', 'GCG', 'GCU'], 
    'GLY': ['GGU', 'GGG', 'GGA', 'GGC'], 
    'ILE': ['AUC', 'AUA', 'AUU'], 
    'LEU': ['UUA', 'UUG', 'CUC', 'CUU', 'CUG', 'CUA'], 
    'HIS': ['CAU', 'CAC'], 
    'ARG': ['CGA', 'CGC', 'CGG', 'CGU', 'AGG', 'AGA'], 
    'TRP': ['UGG'], 
    'VAL': ['GUA', 'GUC', 'GUG', 'GUU'], 
    'GLU': ['GAG', 'GAA'], 
    'TYR': ['UAU', 'UAC']
    }
    
    rscu_values = synonymous_codon_dict.copy()
    
    for aa in synonymous_codon_dict:
        
        codons = synonymous_codon_dict[aa]
        Ni = (len(synonymous_codon_dict[aa]))
        total_count = 0.0
        
        for codon in codons:
            total_count += codon_count_dict[codon]
        
        #calculating RSCU
        
        rscu_values[aa] = {}
        
        for codon in codons:
            rscu = (codon_count_dict[codon] * Ni) / total_count
            rscu_values[aa][codon] = rscu

    
    return rscu_values
            

def count_codon(sequence):
    
    codon_count = codon_dict.copy()
    
    start_pos = useful.get_start(sequence)
    end_pos = useful.get_stop(sequence)
    
    for i in range(start_pos, end_pos-2, 3):
        codon = sequence[i:i+3]
        if codon in codon_count:
            codon_count[codon] += 1
        else:
            raise TypeError ('Unidentified codon {}'.format(codon))       
           
    return codon_count 

def interpret(file, sheet):
    
    count_dict = codon_dict.copy()
    
    df = pd.read_excel(file, sheetname = sheet , index_col = None)
    symbol = df['Symbol'].copy()
    
    ######### find minimum and use as baseline for foldchange #########
    # minimum_intensity = df['Average'].min()
    minimum_intensity = 3.650557279586792
    print (minimum_intensity)
    
    
    for ids in range(len(symbol)):
        print (symbol[ids])
        seq = useful.pull_fasta_sequence(symbol[ids])
        seq = useful.clean_seq(seq)
        seq = Seq.Seq(seq, Seq.Alphabet.generic_dna)
        seq = seq.transcribe()  
        updated_count = count_codon(seq)
        
        fold_change = df['Average'][ids] / minimum_intensity

        updated_count.update((k, v * fold_change) for k, v in updated_count.items())
        count_dict = Counter(count_dict) + Counter(updated_count)
        
    
    rscu_values = get_rscu_value(count_dict)  
    return rscu_values


if __name__ == '__main__':
    print (interpret('tumours_individual_sets_clean.xlsx', 'Ala 20001 - 30000'))