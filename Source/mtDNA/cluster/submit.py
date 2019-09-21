from Source.mtDNA.cluster.rf import random_forest

data_path = 'D:/Aaron/Bio/mtDNA/Data/'
data_path_npz = 'D:/Aaron/Bio/mtDNA/Data/genes/npz/'
data_path_pkl = 'D:/Aaron/Bio/mtDNA/Data/genes/pkl/'
experiment_type = 'mt'
reference_pop = 'IBS'
target_pop = 'FIN'
reference_part = 0.75
result_file_suffix = ''
target_accuracy = 0.68
num_features = 0
gene_files = ['mt_gene_list.txt']
k_mt_max = 1
k_nuc_max = 2

for k_mt in range(1, k_mt_max + 1):
    for k_nuc in range(1, k_nuc_max + 1):
        config_file_name = 'config.txt'
        with open(config_file_name, "w") as f:
            f.write('data_path\t' + data_path + '\n')
            f.write('data_path_npz\t' + data_path_npz + '\n')
            f.write('data_path_pkl\t' + data_path_pkl + '\n')
            f.write('experiment_type\t' + experiment_type + '\n')
            f.write('reference_pop\t' + reference_pop + '\n')
            f.write('target_pop\t' + target_pop + '\n')
            f.write('reference_part\t' + str(reference_part) + '\n')
            f.write('result_file_suffix\t' + result_file_suffix + '\n')
            f.write('target_accuracy\t' + str(target_accuracy) + '\n')
            f.write('num_features\t' + str(num_features) + '\n')
            f.write('gene_files\t' + ', '.join(gene_files) + '\n')
            f.write('k_mt\t' + str(k_mt) + '\n')
            f.write('k_nuc\t' + str(k_nuc) + '\n')
        random_forest()
