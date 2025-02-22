'''
Created on 7 Mar 2018

@author: filipe
''' 

import pandas as pd
import useful

file_out = 'cl_individual_sets_mock.xlsx'
file_name = 'Cell_lines_individual_datasets.xlsx'
sheet = 'Mock'

remove = []
df = pd.read_excel(io = file_name, sheetname = sheet)

#df = df.drop(['Unique ID'], axis = 1)
symbol = df['Symbol'].copy()
df = df.set_index('Symbol')

'''
for every gene in the dataset this will attempt to pull the fasta file
if an error occurs accessing the file, the gene will be added to a remove list
and removed from the dataset

the genes are removed because not all of the gene IDs in the dataset will be endogenous to the species 
as some are addded during the preceding experiments - these will not have fasta sequences on the NCBI database
'''
for i in symbol:
    try:
        useful.pull_fasta_sequence(i)
    except:  
        remove.append(i)
        
df = df.drop(remove)

with pd.ExcelWriter(file_out, engine = 'xlsxwriter') as writer:
    df.to_excel(writer)
    writer.save()
    writer.close()




