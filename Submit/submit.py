import pathlib
import os.path
import itertools
import json
import hashlib
from scipy.special import binom
import sklearn

medium = 1

data_path = '/data/biophys/denysov/yusipov/mtDNA/input/'
data_path_npz = '/data/biophys/denysov/yusipov/mtDNA/input/genes/npz/'
data_path_pkl = '/data/biophys/denysov/yusipov/mtDNA/input/genes/pkl/'
experiment_type = 'mt-nuc'
random_forest_type = 2
reference_pop = 'IBS'
target_pop = 'TSI'
reference_part = 0.75
result_file_suffix = 'diet'
target_accuracy = 0.6
num_features = 1000
gene_files = ['mt_gene_list.txt', 'test_gene_list_diet.txt']
create_tree = 0
run_timer = 0
k_mt_min = 13
k_nuc_min = 9
k_mt_max = 13
k_nuc_max = 9
num_cluster_tasks = 1000
num_atomic_tasks = 50
num_running_tasks = 0

mt_num = 0
nuc_num = 0

for k_mt in range(k_mt_min, k_mt_max + 1):
    mt_num += binom(k_mt_max, k_mt)
for k_nuc in range(k_nuc_min, k_nuc_max + 1):
    nuc_num += binom(k_nuc_max, k_nuc)

num_combinations = mt_num * nuc_num
print('Number of cluster tasks: ' + str(num_cluster_tasks))
print('Number of atomic tasks: ' + str(num_atomic_tasks))
print('Number of combinations: ' + str(int(num_combinations)))

genes_mt = []
genes_nuc = []
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

for k_mt in range(k_mt_min, k_mt_max + 1):
    for k_nuc in range(k_nuc_min, k_nuc_max + 1):
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

if len(result_file_suffix) > 0:
    result_file_suffix = '_' + result_file_suffix

for task_id in range(0, num_cluster_tasks):

    if len(combinations[0]) < (task_id + 1) * num_atomic_tasks:
        genes_mt_task = combinations[0][task_id * num_atomic_tasks:]
        genes_nuc_task = combinations[1][task_id * num_atomic_tasks:]
    else:
        genes_mt_task = combinations[0][task_id * num_atomic_tasks: (task_id + 1) * num_atomic_tasks]
        genes_nuc_task = combinations[1][task_id * num_atomic_tasks: (task_id + 1) * num_atomic_tasks]

    json_list = json.dumps([genes_mt_task, genes_nuc_task]).encode('utf-8')

    curr_hash = hashlib.md5(json_list).hexdigest()

    root = '/data/biophys/denysov/yusipov/mtDNA/output'
    local_path = '/' + experiment_type + '/rf_type_' + str(random_forest_type) + \
                 '/ref_' + reference_pop + '_target_' + target_pop + '/' + \
                 'nat_' + str(num_atomic_tasks) + '/' + curr_hash + '/'
    fn_path = root + local_path
    pathlib.Path(fn_path).mkdir(parents=True, exist_ok=True)

    fn_test = str(target_accuracy) + '_accuracy' + result_file_suffix + '.txt'

    if not os.path.isfile(fn_path + fn_test):

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
        file_config.write('random_forest_type\t' + str(random_forest_type) + '\n')
        file_config.write('reference_pop\t' + reference_pop + '\n')
        file_config.write('target_pop\t' + target_pop + '\n')
        file_config.write('reference_part\t' + str(reference_part) + '\n')
        file_config.write('result_file_suffix\t' + result_file_suffix + '\n')
        file_config.write('target_accuracy\t' + str(target_accuracy) + '\n')
        file_config.write('num_features\t' + str(num_features) + '\n')
        file_config.write('gene_files\t' + ', '.join(gene_files) + '\n')
        file_config.write('create_tree\t' + str(create_tree) + '\n')
        file_config.write('run_timer\t' + str(run_timer) + '\n')

        file_config.close()

        if medium == 0:
            os.system('sbatch run_mpipks_sd_sbatch.sh ' + fn_path)
        elif medium == 1:
            os.system('sbatch run_mpipks_sd_sbatch_medium.sh ' + fn_path)

        num_running_tasks += 1
        if len(combinations[0]) < (task_id + 1) * num_atomic_tasks:
            break

print('Number of running tasks: ' + str(num_running_tasks))
