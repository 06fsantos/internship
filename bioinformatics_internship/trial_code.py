'''
Created on 27 Feb 2018

@author: filipe
'''
import pandas as pd
import useful

def prepare_file(file_name):
    df = pd.read_excel(io = file_name, sheetname = sheet)
    df = df.drop(['Unique ID'], axis = 1)
    symbol = df['Symbol'].copy()
    return symbol



if __name__ == '__main__':
    file = 'tumours.xlsx'
    sheet = 'Mock vs Ala up'
    for i in prepare_file(file):
        try:
            useful.pull_fasta_sequence(i)
        except:
            print ('{} is not a recognised gene_id on NCBI'.format(i))