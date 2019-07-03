'''
Created on 9 Apr 2018

@author: filipe
'''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO, Seq
import useful


def update_dict(file, sheet):
    '''
    Input: the file containing worksheets of the upregulated and downregulated genes, and a dictionary containing 
           all of the codons
    
    output: a dictionary containing the codon frequency per 1000 codons
    '''
    
    
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
        
    codon_count = 0
    
    df = pd.read_excel(file, sheetname = sheet, index_col = None)
    '''    if df.iloc[0]['FC'] > 0:
        df = df.nlargest(n = 10, columns = ['FC'])
    else:
        df = df.nsmallest(n = 10, columns = ['FC'])'''
    symbol = df['Symbol'].copy()
    
    for gene_id in symbol:
        seq = useful.pull_fasta_sequence(gene_id)
        seq = useful.clean_seq(seq)
        seq = Seq.Seq(seq, Seq.Alphabet.generic_dna)
        seq = seq.transcribe()
        start_pos = useful.get_start(seq)
        stop_pos = useful.get_stop(seq)
        
        for j in range(start_pos+3, stop_pos - 2, 3):
            for key in codon_dict:
                if seq[j:j+3] == key:
                    codon_dict[key] += 1
                    codon_count +=1
                    
    print(codon_count)

    codon_dict.update((k, v / codon_count) for k, v in codon_dict.items())
    codon_dict.update((k, round(v,3)) for k, v in codon_dict.items())
    return codon_dict

def plot_codon_frequency(up_dict, down_dict, title):

    names = up_dict.keys()
    values_up = up_dict.values()
    values_down = down_dict.values()
    
    width = 0.4  # width of the bars
    n = len(names)  # the number of codons to be plotted 
    loc = np.arange(n)  # produces evenly spaced values within n
    
    fig = plt.figure(figsize = (15.0, 7.0))
    ax = fig.add_subplot(111)
    
    up = ax.bar(loc, values_up, width = width, color = 'r', alpha = 1.0)
    down = ax.bar(loc + width, values_down, width = width, color = 'k', alpha = 1.0)
    
    ax.set_title(title)
    ax.set_ylabel('Codon Frequency(/Total number of codons)')
    ax.set_xlabel('Codon')
    ax.set_xticks(loc + width / 2)
    ax.set_xticklabels(names, rotation = 45, rotation_mode = 'anchor', ha = 'right')
    ax.legend((up[0], down[0]), ('Up', 'Down'))
    
    plt.show()

if __name__ == '__main__':
    file = 'tumours_clean.xlsx'
    title = 'Tumours: Mock vs Wt'
    sheet_up = 'Mock vs Wt up'
    sheet_down = 'Mock vs Wt down'
    plot_codon_frequency(update_dict(file, sheet_up), update_dict(file, sheet_down), title)