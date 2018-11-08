import os
import numpy as np

data_path = '../Data/'
count_path = data_path + 'snp_counts/'

gene_data = []
snp_data = []
chromosome_data = []
for chr_data_file in os.listdir(count_path):
    chr_name = chr_data_file[:-14]
    f = open(count_path + chr_data_file)
    for line in f:
        line = line.replace('\n', '')
        curr_gene_data = line.split(' ')
        gene_name = curr_gene_data[0]
        snp_count = int(curr_gene_data[1])
        gene_data.append(gene_name)
        snp_data.append(snp_count)
        chromosome_data.append(chr_name)
    f.close()

order = np.argsort(snp_data)
snp_data = list(np.array(snp_data)[order])
gene_data = list(np.array(gene_data)[order])
chromosome_data = list(np.array(chromosome_data)[order])
fn = data_path + 'gene_snp_statistics.txt'
data = list(map(list, zip(*[gene_data, chromosome_data, snp_data])))
np.savetxt(fn, data, fmt='%s')
