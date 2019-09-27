import pathlib
import os.path
import itertools
import json
import hashlib
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

for file_id in range(0, len(gene_files)):
    data_gene_file = open(data_path + gene_files[file_id])
    if file_id == 0:
        if experiment_type == 'mt':
            genes_mt = [i for i, line in enumerate(data_gene_file)]
            genes_nuc = [[]]
        elif experiment_type == 'nuc':
            genes_mt = [[]]
            genes_nuc = [i for i, line in enumerate(data_gene_file)]
        elif experiment_type == 'mt-nuc':
            genes_mt = [i for i, line in enumerate(data_gene_file)]
    else:
        genes_nuc = [i for i, line in enumerate(data_gene_file)]
    data_gene_file.close()

combinations = [[], []]

for k_mt in range(1, k_mt_max + 1):
    for k_nuc in range(1, k_nuc_max + 1):
        print('k_mt: ' + str(k_mt))
        print('k_nuc: ' + str(k_nuc))

        for subset_mt in itertools.combinations(genes_mt, k_mt):
            for subset_nuc in itertools.combinations(genes_nuc, k_nuc):
                if experiment_type == 'mt':
                    combinations[0].append(list(subset_mt))
                    combinations[1].append(list(subset_nuc)[0])
                elif experiment_type == 'nuc':
                    combinations[0].append(list(subset_mt)[0])
                    combinations[1].append(list(subset_nuc))
                else:
                    combinations[0].append(list(subset_mt))
                    combinations[1].append(list(subset_nuc))

for task_id in range(0, num_cluster_tasks):

    if len(combinations[0]) < (task_id + 1) * num_atomic_tasks:
        genes_mt_task = combinations[0][task_id * num_atomic_tasks :]
        genes_nuc_task = combinations[1][task_id * num_atomic_tasks :]
    else:
        genes_mt_task = combinations[0][task_id * num_atomic_tasks : (task_id + 1) * num_atomic_tasks]
        genes_nuc_task = combinations[1][task_id * num_atomic_tasks : (task_id + 1) * num_atomic_tasks]

    json_list = json.dumps([genes_mt_task, genes_nuc_task]).encode('utf-8')

    hash = hashlib.md5(json_list).hexdigest()

    root = '/data/biophys/denysov/yusipov/mtDNA/output'
    local_path = '/' + experiment_type + '/ref_' + reference_pop + '_target_' + target_pop + '/' + hash + '/'
    fn_path = root + local_path
    pathlib.Path(fn_path).mkdir(parents=True, exist_ok=True)

    file_mt = open(fn_path + '/config_mt_genes.txt', 'w')
    for genes_task in genes_mt_task:
        file_mt.write('\t'.join([str(item) for item in genes_task]) + '\n')
    file_mt.close()

    file_nuc = open(fn_path + '/config_nuc_genes.txt', 'w')
    for genes_task in genes_nuc_task:
        file_nuc.write('\t'.join([str(item) for item in genes_task]) + '\n')
    file_nuc.close()

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

    file_config.close()

    if medium == 0:
        os.system('sbatch run_mpipks_sd_sbatch.sh ' + fn_path)
    elif medium == 1:
        os.system('sbatch run_mpipks_sd_sbatch_medium.sh ' + fn_path)
