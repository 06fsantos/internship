'''
Created on 28 Feb 2018

@author: filipe
'''
from Bio import Entrez, SeqIO, Seq


def pull_fasta_sequence(gene_id):
    '''
    pull the gene sequence from the NCBI nucleotide databse
    
    -------------------------------
    input:
        gene_id:
            the ID of the gene whose sequence is being pulled from the NCBI database 
        
    returns:
        the nucleotide sequence of the gene 
    '''
    Entrez.email = 'santos.filipe1995@gmail.com'
    with Entrez.efetch(db = 'nucleotide', id = gene_id, rettype = 'fasta', retmode = 'text') as handle:
        seq_record = SeqIO.read(handle, 'fasta')
        return (seq_record.seq)

def clean_seq(sequence):
    '''
    removes any letters unless ACTG from the sequence ensures no trailing information 
    attached to the end of the sequence is misinterpreted as a part of the sequence 
    
    --------------------------------
    input:
        sequence: the DNA sequence to be cleaned 
            
    returns:
        seq_return: the cleaned sequence with only the bases A, T, G and C
    '''
    DNA_BASES = ['A', 'T', 'C', 'G']
    seq_return = ''
    for i in sequence:
        if i.upper() in DNA_BASES:
            seq_return += i.upper()
    return seq_return 
   
def get_start(sequence):
    '''
    Identifies the start codon given an RNA sequence as the input 
    
    -----------------------------
    Input:
        sequence: the RNA sequence to be analysed
        
    Returns:
        start: the location of the first start codon within the sequence 
    '''
    start_codon = 'AUG'
    for i in range(0, len(sequence)):
        if sequence[i:i+3] == start_codon:
            start = i
            break
    return start

def get_stop(sequence):
    '''
    Identifies the position of the stop codon 
    
    --------------------------------
    Input:
        sequence: the RNA sequence to be analysed
    
    Returns:
        stop: the location of the first stop codon that follows the start codon in the sequence
    '''
    stop_codons = ['UAA', 'UAG', 'UGA']
    for i in range(get_start(sequence), len(sequence) -2, 3):
        if sequence[i:i+3] in stop_codons:
            stop = i
            break
        else:
            stop = len(sequence)
    return stop
    
def codon_percentage(sequence, codon):
    '''
    counts the number of the desired codon and the total number of codons and calculates a percentage from both
    does not include START or STOP codons in the total_codon count 
    
    --------------------------------------
    input: 
        sequence: the desired sequence to analyse
    
    Returns: 
        percentage: the percentage of codons within the sequence that are the desired codon
    
    '''
    stop_codon = ['UAA', 'UAG', 'UGA']
    codon_count = 0
    total_codon = 0
    seq = clean_seq(sequence)
    seq = Seq.Seq(seq, Seq.Alphabet.generic_dna)
    seq = seq.transcribe()
    start_pos = get_start(seq)
    for j in range(start_pos+3, len(seq) - 2, 3):
        if seq[j:j+3] == codon:
            codon_count += 1
            total_codon += 1
        else:
            total_codon += 1
            if seq[j:j+3] in stop_codon:
                total_codon -= 1
                break
    percentage = (codon_count/total_codon) * 100
    return percentage