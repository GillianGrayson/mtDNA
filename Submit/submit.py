import pathlib
import os.path
import numpy as np
import sklearn

medium = 0

data_path = '/data/biophys/denysov/yusipov/mtDNA/input/'
data_path_npz = '/data/biophys/denysov/yusipov/mtDNA/input/genes/npz/'
data_path_pkl = '/data/biophys/denysov/yusipov/mtDNA/input/genes/pkl/'
experiment_type = 'mt'
reference_pop = 'IBS'
target_pop = 'FIN'
reference_part = 0.75
result_file_suffix = ''
target_accuracy = 0.68
num_features = 0
gene_files = ['mt_gene_list.txt']
k_mt_max = 1
k_nuc_max = 1

for k_mt in range(1, k_mt_max + 1):
    for k_nuc in range(1, k_nuc_max + 1):

        print('k_mt: ' + str(k_mt))
        print('k_nuc: ' + str(k_nuc))

        root = '/data/biophys/denysov/yusipov/mtDNA/output'
        local_path = '/' + experiment_type  + '/ref(' + reference_pop + ')_target(' + target_pop + ')/' + 'mt(' + str(k_mt) + ')_nuc(' + str(k_nuc) + ')'
        fn_path = root + local_path
        pathlib.Path(fn_path).mkdir(parents=True, exist_ok=True)

        file_config = open(fn_path + '/config.txt', 'w')

        file_config.write('data_path\t' + data_path + '\n')
        file_config.write('data_path_npz\t' + data_path_npz + '\n')
        file_config.write('data_path_pkl\t' + data_path_pkl + '\n')
        file_config.write('experiment_type\t' + experiment_type + '\n')
        file_config.write('reference_pop\t' + reference_pop + '\n')
        file_config.write('target_pop\t' + target_pop + '\n')
        file_config.write('reference_part\t' + str(reference_part) + '\n')
        file_config.write('result_file_suffix\t' + result_file_suffix + '\n')
        file_config.write('target_accuracy\t' + str(target_accuracy) + '\n')
        file_config.write('num_features\t' + str(num_features) + '\n')
        file_config.write('gene_files\t' + ', '.join(gene_files) + '\n')
        file_config.write('k_mt\t' + str(k_mt) + '\n')
        file_config.write('k_nuc\t' + str(k_nuc) + '\n')

        file_config.close()

        # if medium == 0:
        #     os.system('sbatch run_mpipks_sd_sbatch.sh ' + fn_path)
        # elif medium == 1:
        #     os.system('sbatch run_mpipks_sd_sbatch_medium.sh ' + fn_path)
        
        print('Hello')
