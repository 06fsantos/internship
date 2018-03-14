'''
Created on 9 Mar 2018

@author: filipe
'''
import pandas as pd
import numpy as np
import useful
import matplotlib.pyplot as plt

plt.style.use('ggplot')

def percentages(file, sheet, codon):
    df = pd.read_excel(file, sheetname = sheet, index_col = None)
    symbol = df['Symbol'].copy()
    for codon in codons:
        codon_percent = []
        for gene_id in symbol:
            seq = useful.pull_fasta_sequence(gene_id)
            percentage = useful.codon_percentage(seq, codon)
            codon_percent.append(percentage)
        df[codon] = codon_percent
    print (df)
    return df

def plot_avg_percent(dataframe_up, dataframe_down, title):
    '''
    input: dataframe containing codon percentage in the upregulated genes and the downregulated genes,
           also the title of the output plot
    
    output: a bar chart displaying differences in codon bias within with up and downregulated genes
    '''
    mean_dict_up = {}
    sem_dict_up = {}
    for column in dataframe_up:
        if column in codons:
            mean = dataframe_up[column].mean()
            std_err = dataframe_up[column].sem()
            mean_dict_up[column] = mean
            sem_dict_up[column] = std_err
    
    mean_dict_down = {}
    sem_dict_down = {}
    for column in dataframe_down:
        if column in codons:
            mean = dataframe_down[column].mean()
            std_err = dataframe_down[column].sem()
            mean_dict_down[column] = mean
            sem_dict_down[column] = std_err
    
    names = mean_dict_up.keys()
    values_up = mean_dict_up.values()
    errors_up = sem_dict_up.values()
    
    values_down = mean_dict_down.values()
    errors_down = sem_dict_down.values()

    width = 0.35 # width of the bars
    n = len(names) # number of distinct desired codons to be plotted
    loc = np.arange(n)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    up = ax.bar(loc, values_up, width = width, color = 'g', yerr = errors_up, alpha = 0.6)
    down = ax.bar(loc + width, values_down, width = width, color = 'b', yerr = errors_down, alpha = 0.6)
    
    ax.set_title(title)
    ax.set_ylabel('codon percentages')
    ax.set_xlabel('codons')
    ax.set_xticks(loc + width / 2)
    ax.set_xticklabels(names)
    ax.legend((up[0], down[0]), ('up', 'down'))
    
    plt.show()
        

if __name__ == '__main__':
    file = 'tumours_clean.xlsx'
    title = 'Tumours: Mock vs Ala'
    sheet_up = 'Mock vs Ala up'
    sheet_down = 'Mock vs Ala down'
    codons = ['GCU', 'GCC', 'GCA']
    plot_avg_percent(percentages(file, sheet_up, codons), percentages(file, sheet_down, codons), title)