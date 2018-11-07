import os
import numpy as np

data_path = '../Data/'
target_chromosome = 'chr20'
result_path = '../Data/' + target_chromosome + '/'

gene_data = []
for gene_name in os.listdir(result_path):
    f = open(result_path + gene_name)
    num_snps = 0
    for line in f:
        num_snps += 1
    gene_data.append([gene_name[:-4], num_snps-1])
fn = data_path + target_chromosome + '_snp_count.txt'
np.savetxt(fn, gene_data, fmt='%s')
