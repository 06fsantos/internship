'''
Created on 11 Apr 2018

@author: filipe
'''
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO, Seq
import useful

def update_dict(gene_id):
    '''
    isolates every codon in a DNA sequence and returns a dictionary containing the count of each codon
    
    --------------------------------
    Input: 
        gene_id: 
            The Gene_id of the target sequence
    
    Returns: 
        codon_dict:
            a dictionary containing the codon composition of the gene
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


    return codon_dict

def plot_codon_composition(composition_dict, title):
    '''
    produces a plot of the codon composition of a DNA sequence 
    
    ------------------------------------
    Input:
        composition_dict:
            a dictionary containing the codon frequency for every codon 
        
        title:
            the title of the plot
    
    Output:
        produces and displays a bar chart denoting the codon composition of the DNA sequence 
    '''
    names = composition_dict.keys()
    values = composition_dict.values()
    
    width = 0.8
    n = len(names)
    loc = np.arange(n)
    
    fig = plt.figure(figsize = (15.0, 7.0))
    ax = fig.add_subplot(111)
    
    ax.bar(loc, values, width = width, color = 'r')
    ax.set_title(title)
    ax.set_ylabel('Number of Codons')
    ax.set_xlabel('Codon')
    ax.set_xticks(loc)
    ax.set_xticklabels(names, rotation = 45, rotation_mode = 'anchor', ha = 'right')
    
    plt.show()

if __name__ == '__main__':
    gene = ['XM_892296', 'NM_001162901', 'NM_001033330']
    for i in gene:
        title = 'Codon Compostion: {} (Cell Lines Mock vs Wt down)'.format(i)
        plot_codon_composition(update_dict(i), title)