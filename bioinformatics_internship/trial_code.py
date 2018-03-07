'''
Created on 27 Feb 2018

@author: filipe
'''
import pandas as pd
import useful

def prepare_file(file_name):
    remove = []
    df = pd.read_excel(io = file_name, sheetname = sheet)
    df = df.drop(['Unique ID'], axis = 1)
    symbol = df['Symbol'].copy()
    df = df.set_index('Symbol')
    for i in symbol:
        try:
            useful.pull_fasta_sequence(i)
        except: 
            remove.append(i)
    df = df.drop(remove)
    return df


if __name__ == '__main__':
    file = 'tumours.xlsx'
    sheet = 'Mock vs Ala up'
    print (prepare_file(file))


        