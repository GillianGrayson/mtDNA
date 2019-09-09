import os
import numpy as np

data_path = '../Data/'
target_chromosome = 'chrX'
result_path = data_path + target_chromosome + '/'
count_path = data_path + '/snp_counts/'
if not os.path.exists(count_path):
    os.makedirs(count_path)

gene_data = []
for gene_file_name in os.listdir(result_path):
    gene_name = gene_file_name[:-4]
    f = open(result_path + gene_file_name)
    num_snps = 0
    for line in f:
        num_snps += 1
    gene_data.append([gene_name, num_snps-1])
    f.close()
fn = count_path + target_chromosome + '_snp_count.txt'
np.savetxt(fn, gene_data, fmt='%s')
