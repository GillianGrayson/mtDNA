import json
import hashlib
import pathlib
import os

medium = 1

data_path = '/data/biophys/denysov/yusipov/mtDNA/input/'
data_path_npz = '/data/biophys/denysov/yusipov/mtDNA/input/genes/npz/'
data_path_pkl = '/data/biophys/denysov/yusipov/mtDNA/input/genes/pkl/'
experiment_type = 'mt-nuc'
random_forest_type = 2
reference_pop = 'FIN'
target_pop = 'IBS'
reference_part = 0.75
result_file_suffix = ''
target_accuracy = 0.6
num_features = 0
gene_files = ['test_mt.txt', 'test_nuc.txt']
create_tree = 0
run_timer = 0
num_cluster_tasks = 1
num_atomic_tasks = 100
num_running_tasks = 0

result_path = '/data/biophys/denysov/yusipov/mtDNA/output/'
experiment_result_path = result_path + experiment_type + '/' + \
                         'ref_' + reference_pop + '_' + 'target_' + target_pop + '/'

genes_mt = []
genes_mt_names = []
genes_nuc = []
genes_nuc_names = []
for file_id in range(0, len(gene_files)):
    data_gene_file = open(data_path + gene_files[file_id])
    if file_id == 0:
        if experiment_type == 'mt':
            for i, line in enumerate(data_gene_file):
                genes_mt.append(i)
                genes_mt_names.append(line.rstrip())
            genes_nuc = []
        elif experiment_type == 'nuc':
            genes_mt = []
            for i, line in enumerate(data_gene_file):
                genes_nuc.append(i)
                genes_nuc_names.append(line.rstrip())
        elif experiment_type == 'mt-nuc':
            for i, line in enumerate(data_gene_file):
                genes_mt.append(i)
                genes_mt_names.append(line.rstrip())
    else:
        for i, line in enumerate(data_gene_file):
            genes_nuc.append(i)
            genes_nuc_names.append(line.rstrip())
    data_gene_file.close()

if len(result_file_suffix) > 0:
    result_file_suffix = '_' + result_file_suffix

json_list = json.dumps([[genes_mt], [genes_nuc]]).encode('utf-8')

curr_hash = hashlib.md5(json_list).hexdigest()

root = result_path
local_path = '/' + experiment_type + '/rf_type_' + str(
    random_forest_type) + '/ref_' + reference_pop + '_target_' + target_pop + '/' + curr_hash + '/'
fn_path = root + local_path
pathlib.Path(fn_path).mkdir(parents=True, exist_ok=True)

experiment_type_suffix = 'top_features_'
if experiment_type == 'mt':
    experiment_type_suffix += 'mt_' + '_'.join(list(map(str, genes_mt))) + result_file_suffix
elif experiment_type == 'nuc':
    experiment_type_suffix += 'nuc_' + '_'.join(list(map(str, genes_nuc))) + result_file_suffix
elif experiment_type == 'mt-nuc':
    experiment_type_suffix += 'mt_' + '_'.join(list(map(str, genes_mt))) + \
                              '_nuc_' + '_'.join(list(map(str, genes_nuc))) + result_file_suffix

fn_features = str(target_accuracy) + '_' + experiment_type_suffix + result_file_suffix + '.txt'
top_features = []

features_file = open(fn_path + fn_features)
for line in features_file:
    top_features.append(line.split('\t')[0])
features_file.close()

combinations = [[], []]
if experiment_type == 'mt':
    top_features = [genes_mt_names.index(top_features[i]) for i in range(0, len(top_features))]
    combinations[0] = [top_features[:i] for i in range(1, len(top_features))]
    combinations[1].append([genes_nuc])
elif experiment_type == 'nuc':
    top_features = [genes_nuc_names.index(top_features[i]) for i in range(0, len(top_features))]
    combinations[0].append([genes_mt])
    combinations[1] = [top_features[:i] for i in range(1, len(top_features))]
else:
    mt_genes = []
    nuc_genes = []
    for i in range(0, len(top_features)):
        curr_features = top_features[i].split('_')
        if len(curr_features) == 2:
            curr_mt_gene = genes_mt_names.index(curr_features[0])
            curr_nuc_gene = genes_nuc_names.index(curr_features[1])
        else:
            curr_mt_gene = genes_mt_names.index(curr_features[0] + '_' + curr_features[1])
            curr_nuc_gene = genes_nuc_names.index(curr_features[2])
        mt_genes.append(curr_mt_gene)
        nuc_genes.append(curr_nuc_gene)
        combinations[0].append(sorted(set(mt_genes[:i+1]), key=mt_genes[:i+1].index))
        combinations[1].append(sorted(set(nuc_genes[:i+1]), key=nuc_genes[:i+1].index))

for task_id in range(0, num_cluster_tasks):

    if len(combinations[0]) < (task_id + 1) * num_atomic_tasks:
        genes_mt_task = combinations[0][task_id * num_atomic_tasks:]
        genes_nuc_task = combinations[1][task_id * num_atomic_tasks:]
    else:
        genes_mt_task = combinations[0][task_id * num_atomic_tasks: (task_id + 1) * num_atomic_tasks]
        genes_nuc_task = combinations[1][task_id * num_atomic_tasks: (task_id + 1) * num_atomic_tasks]

    json_list = json.dumps([genes_mt_task, genes_nuc_task]).encode('utf-8')

    curr_hash = hashlib.md5(json_list).hexdigest()

    root = 'D:/Aaron/Bio/mtDNA/Result/files'
    local_path = '/' + experiment_type + '/rf_type_' + str(random_forest_type) + \
                 '/ref_' + reference_pop + '_target_' + target_pop + '/' + curr_hash + '/'
    fn_path = root + local_path
    pathlib.Path(fn_path).mkdir(parents=True, exist_ok=True)

    fn_test = str(target_accuracy) + '_accuracy.txt'

    if not os.path.isfile(fn_test):

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
