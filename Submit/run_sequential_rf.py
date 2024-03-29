import json
import hashlib
import pathlib
import os
import time


medium = 1

data_path = '/data/condmat/ivanchen/kalyakulina/mtDNA/input/'
data_path_npz = '/data/condmat/ivanchen/kalyakulina/mtDNA/input/genes/npz/'
data_path_pkl = '/data/condmat/ivanchen/kalyakulina/mtDNA/input/genes/pkl/'
experiment_type = 'mt'
random_forest_type = 3
reference_pop = 'GBR'
target_pop = 'FIN'
reference_part = 0.75
result_file_suffix = ''
target_accuracy = 0.55
num_features = 0
num_sequential_runs = 1
sequential_run_type = 'lin'  # lin, max, exp
num_estimators = 500
num_cv_runs = 5
gene_files = ['mt_gene_list.txt']
create_tree = 0
run_timer = 1
num_cluster_tasks = 1
num_atomic_tasks = 1

result_path = '/data/condmat/ivanchen/kalyakulina/mtDNA/output/'

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

json_list = json.dumps([genes_mt, genes_nuc]).encode('utf-8')

curr_hash = hashlib.md5(json_list).hexdigest()

root = '/data/condmat/ivanchen/kalyakulina/mtDNA/output/'
local_path = '/' + experiment_type + '/rf_type_' + str(random_forest_type) + \
             '/ref_' + reference_pop + '_target_' + target_pop + '/' + \
             'nat_' + str(num_atomic_tasks) + '/' + curr_hash + '/'
fn_path = root + local_path
pathlib.Path(fn_path).mkdir(parents=True, exist_ok=True)

fn_test = str(target_accuracy) + '_accuracy.txt'

if not os.path.isfile(fn_test):

    file_mt = open(fn_path + '/config_mt_genes.txt', 'w')
    file_mt.write('\t'.join([str(item) for item in genes_mt]) + '\n')
    file_mt.close()

    file_nuc = open(fn_path + '/config_nuc_genes.txt', 'w')
    file_nuc.write('\t'.join([str(item) for item in genes_nuc]) + '\n')
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
    file_config.write('num_sequential_runs\t' + str(num_sequential_runs) + '\n')
    file_config.write('sequential_run_type\t' + sequential_run_type + '\n')
    file_config.write('num_estimators\t' + str(num_estimators) + '\n')
    file_config.write('num_cv_runs\t' + str(num_cv_runs) + '\n')
    file_config.write('gene_files\t' + ', '.join(gene_files) + '\n')
    file_config.write('create_tree\t' + str(create_tree) + '\n')
    file_config.write('run_timer\t' + str(run_timer) + '\n')
    file_config.close()

    if run_timer == 1:
        start = time.process_time()

    if medium == 0:
        os.system('sbatch run_mpipks_sd_sbatch.sh ' + fn_path)
    elif medium == 1:
        os.system('sbatch run_mpipks_sd_sbatch_medium.sh ' + fn_path)
