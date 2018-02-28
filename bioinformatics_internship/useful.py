'''
Created on 28 Feb 2018

@author: filipe
'''
from Bio import Entrez, SeqIO, Seq


def pull_fasta_sequence(gene_id):
    '''
    pull the gene sequence from the NCBI nucleotide databse
    '''
    Entrez.email = 'santos.filipe1995@gmail.com'
    with Entrez.efetch(db = 'nucleotide', id = gene_id, rettype = 'fasta', retmode = 'text') as handle:
        seq_record = SeqIO.read(handle, 'fasta')
        return (seq_record.seq)

def clean_seq(sequence):
    '''
    removes any letters unless ACTG
    '''
    DNA_BASES = ['A', 'T', 'C', 'G']
    seq_return = ''
    for i in sequence:
        if i.upper() in DNA_BASES:
            seq_return += i.upper()
    return seq_return 
   

def count_codons(sequence, codon):
    start_codon = 'AUG'
    stop_codon = ['UAA', 'UAG', 'UGA']
    codon_count = 0
    seq = clean_seq(sequence)
    seq = Seq.Seq(seq, Seq.Alphabet.generic_dna)
    seq = seq.transcribe()
    for i in range(0, len(seq)):
        if seq[i:i+3] == start_codon:
            start_pos = i
            for j in range(start_pos, len(seq) - 2, 3):
                if seq[j:j+3] in stop_codon:
                    break
                elif seq[j:j+3] == codon:
                    codon_count += 1
    return codon_count 


                
        
        
        
        
        
        
        
        
