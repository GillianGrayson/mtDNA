import pathlib
import os.path
import itertools
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
target_accuracy = 0.55
num_features = 0
gene_files = ['mt_gene_list.txt']
create_tree = 0
k_mt_max = 1
k_nuc_max = 1
num_cluster_tasks = 100
num_atomic_tasks = 100
genes_mt = []
genes_nuc = []

for file_id in range(0, len(gene_files)):
    data_gene_file = open(data_path + gene_files[file_id])
    if file_id == 0:
        if experiment_type == 'mt':
            genes_mt.append([i for i, line in enumerate(data_gene_file)])
            genes_nuc.append([])
        elif experiment_type == 'nuc':
            genes_mt.append([])
            genes_nuc.append([i for i, line in enumerate(data_gene_file)])
        elif experiment_type == 'mt-nuc':
            genes_mt.append([i for i, line in enumerate(data_gene_file)])
    else:
        genes_nuc.append([i for i, line in enumerate(data_gene_file)])
    data_gene_file.close()

for k_mt in range(1, k_mt_max + 1):
    for k_nuc in range(1, k_nuc_max + 1):
        print('k_mt: ' + str(k_mt))
        print('k_nuc: ' + str(k_nuc))

        for subset_mt in itertools.combinations(genes_mt, k_mt):
            for subset_nuc in itertools.combinations(genes_nuc, k_nuc):

        root = '/data/biophys/denysov/yusipov/mtDNA/output'
        local_path = '/' + experiment_type + '/ref_' + reference_pop + '_target_' + target_pop + '/' + 'mt_' + str(k_mt) + '_nuc_' + str(k_nuc)
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
        file_config.write('create_tree\t' + str(create_tree) + '\n')
        file_config.write('k_mt\t' + str(k_mt) + '\n')
        file_config.write('k_nuc\t' + str(k_nuc) + '\n')

        file_config.close()

        if medium == 0:
            os.system('sbatch run_mpipks_sd_sbatch.sh ' + fn_path)
        elif medium == 1:
            os.system('sbatch run_mpipks_sd_sbatch_medium.sh ' + fn_path)
