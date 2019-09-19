from Source.mtDNA.cluster.rf import random_forest

data_path = '../../../Data/'
experiment_type = 'mt'
reference_pop = 'FIN'
target_pop = 'GBR'
reference_part = 0.75
result_file_suffix = 'mt'
target_accuracy = 0.8
num_features = 0
gene_files = ['mt_gene_list.txt']
k_mt_max = 13
k_nuc_max = 1

config_file_name = 'config.txt'
with open(config_file_name, "w") as f:
    f.write('experiment_type\t' + experiment_type + '\n')
    f.write('reference_pop\t' + reference_pop + '\n')
    f.write('target_pop\t' + target_pop + '\n')
    f.write('reference_part\t' + str(reference_part) + '\n')
    f.write('result_file_suffix\t' + result_file_suffix + '\n')
    f.write('target_accuracy\t' + str(target_accuracy) + '\n')
    f.write('num_features\t' + str(num_features) + '\n')
    f.write('gene_files\t' + ', '.join(gene_files) + '\n')

for k_mt in range(1, k_mt_max + 1):
    for k_nuc in range(1, k_nuc_max + 1):
        random_forest(k_mt, k_nuc)
