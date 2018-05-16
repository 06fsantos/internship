'''
Created on 7 Mar 2018

@author: filipe
''' 
import pandas as pd
import useful

file_out = 'full_genome_clean.xlsx'
file_name = 'full_genome.xlsx'



remove = []
df = pd.read_excel(io = file_name)
print (df)
#df = df.drop(['Unique ID'], axis = 1)
symbol = df['Symbol'].copy()
df = df.set_index('Symbol')
for i in symbol:
    try:
        useful.pull_fasta_sequence(i)
    except:  
        print (i)
        remove.append(i)
df = df.drop(remove)
with pd.ExcelWriter(file_out, engine = 'xlsxwriter') as writer:
    df.to_excel(writer)
    writer.save()
    writer.close()


   

