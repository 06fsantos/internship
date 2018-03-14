'''
Created on 27 Feb 2018

@author: filipe
'''
import pandas as pd
import numpy as np
import useful
import matplotlib.pyplot as plt 

plt.style.use('ggplot')


def percentages(file, sheet, codon):
    df = pd.read_excel(file, sheetname = sheet, index_col = None)
    if df.iloc[0]['FC'] > 0:
        df = df.nlargest(n = int((len(df)+1) * 0.1), columns = ['FC'])
    else:
        df = df.nsmallest(n = int((len(df)+1) * 0.1), columns = ['FC'])
    symbol = df['Symbol'].copy()
    codon_percent = []
    for gene_id in symbol:
        seq = useful.pull_fasta_sequence(gene_id)
        percentage = useful.codon_percentage(seq, codon)
        codon_percent.append(percentage)
    df[codon] = codon_percent
    print (df.describe())
    return df


'''def plot_avg_percent(dataframe_up, dataframe_down):
    '''
'''    input: dataframe containing codon percentage in the upregulated genes and the downregulated genes
    
    output: a bar chart displaying differences in codon bias within with up and downregulated genes'''
'''
    mean_dict_up = {}
    sem_dict_up = {}
    for column in dataframe_up:
        mean = dataframe_up[column].mean()
        std_err = dataframe_up[column].sem()
        mean_dict_up[column] = mean
        sem_dict_up[column] = std_err
    
    mean_dict_down = {}
    sem_dict_down = {}
    for column in dataframe_down:
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
    
    ax.set_title('trial codon histogram')
    ax.set_ylabel('codon percentages')
    ax.set_xlabel('codons')
    ax.set_xticks(loc + width / 2)
    ax.set_xticklabels(names)
    ax.legend((up[0], down[0]), ('up', 'down'))
    
    plt.show()'''
        

if __name__ == '__main__':
    GCU_up = [1.15, 1.47, 2.54, 2.03, 1.76, 2.99, 1.63, 2.83, 1.47]
    GCU_down = [0.23, 0.45, 0.56, 1.02, 1.14, 0.11, 0.89, 1.34, 0.67]
    GCC_up = [1.23, 1.45, 3.46, 4.21, 1.65, 0.45, 1.23, 1.98, 1.72]
    GCC_down = [0.24, 0.74, 1.56, 1.74, 0.68, 0.33, 0.25, 1.81, 1.11]
    GCA_up = [1.39, 0.50, 2.89, 1.06, 1.83, 3.64, 1.07, 1.58, 1.19]
    GCA_down = [0.78, 0.37, 1.25, 1.09, 0.92, 0.22, 1.89, 0.39, 1.03]
    
    df_up = pd.DataFrame({'GCU': GCU_up, 'GCC': GCC_up, 'GCA': GCA_up})
    df_down = pd.DataFrame({'GCU': GCU_down, 'GCC': GCC_down, 'GCA': GCA_down})
    
    file = 'tumours_clean.xlsx'
    title = 'Tumour Cells: Mock vs Ala'
    sheet_up = 'Mock vs Wt up'
    sheet_down = 'Mock vs Wt down'
    codons = ['GCU', 'GCC', 'GCA']
    codon = 'GCU'
    print (percentages(file, sheet_up, codon))
        