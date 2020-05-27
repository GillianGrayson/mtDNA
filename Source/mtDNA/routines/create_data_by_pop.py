import os
import pandas as pd

data_path = 'E:/YandexDisk/mtDNA/Data/'
genes_path = data_path + 'genes/txt/'
target_pops = ['GBR', 'FIN', 'TSI']
gene_file = 'test_gene_list_cold_adaptation.txt'
pop_file = 's_pop.txt'

result_path = data_path + 'cold_adaptation/'
if not os.path.exists(result_path):
    os.makedirs(result_path)

genes = []
f = open(data_path + gene_file, 'r')
for line in f:
    genes.append(line.rstrip())
f.close()

subjects_to_remove = []
subjects_info = {'subject': [], 'pop': []}
f = open(data_path + pop_file)
f.readline()
for line in f:
    line_list = line.rstrip().split('\t')
    if line_list[1] not in target_pops:
        subjects_to_remove.append(line_list[0])
    else:
        subjects_info['subject'].append(line_list[0])
        subjects_info['pop'].append(line_list[1])
f.close()

subjects_df = pd.DataFrame(subjects_info)
subjects_df.to_csv(result_path + 'populations.txt', index=None, sep='\t', mode='a')

for gene in genes:
    data = pd.read_csv(genes_path + gene + '.txt', sep=' ')
    data = data.drop(columns=subjects_to_remove)
    data.to_csv(result_path + gene + '.txt', index=None, sep='\t', mode='a')
