from Source.mtDNA.cluster.random_forest_local import random_forest
import json
import hashlib
import pathlib
import os

data_path = 'D:/Aaron/Bio/mtDNA/Data/'
data_path_npz = 'D:/Aaron/Bio/mtDNA/Data/genes/npz/'
data_path_pkl = 'D:/Aaron/Bio/mtDNA/Data/genes/pkl/'
experiment_type = 'nuc'
random_forest_type = 1
reference_pop = 'IBS'
target_pop = 'FIN'
reference_part = 0.75
result_file_suffix = 'short'
target_accuracy = 0.6
num_features = 100
gene_files = ['test_gene_list_short.txt']
create_tree = 0
run_timer = 0
k_mt_max = 1
k_nuc_max = 1
num_cluster_tasks = 1
num_atomic_tasks = 10
num_running_tasks = 0
num_top_results = 5

result_path = 'D:/Aaron/Bio/mtDNA/Result/cluster/'
experiment_result_path = result_path + experiment_type + '/' + \
                         'ref_' + reference_pop + '_' + 'target_' + target_pop + '/'

top_accuracy = []
top_indices = []
top_dirs = []
top_genes_mt = []
top_genes_nuc = []

min_accuracy_value = 0.0
min_accuracy_index = num_top_results

dirnames = [x[0] for x in os.walk(experiment_result_path)][1:]

for dirname in dirnames:
    for filename in os.listdir(dirname):
        if 'accuracy' in filename:
            f = open(dirname + '/' + filename)
            line_num = 0
            for line in f:
                item = line.replace('\n', '')
                if len(top_accuracy) < num_top_results:
                    top_accuracy.append(float(item))
                    min_accuracy_value = min(top_accuracy)
                    min_accuracy_index = top_accuracy.index(min_accuracy_value)
                    top_indices.append(line_num)
                    top_dirs.append(dirname)
                else:
                    if float(item) > min_accuracy_value:
                        top_accuracy[min_accuracy_index] = float(item)
                        min_accuracy_value = min(top_accuracy)
                        min_accuracy_index = top_accuracy.index(min_accuracy_value)
                        top_indices[min_accuracy_index] = line_num
                        top_dirs[min_accuracy_index] = dirname
                line_num += 1
            f.close()

if experiment_type == 'mt':
    for id in range(0, len(top_dirs)):
        for file in os.listdir(top_dirs[id]):
            if 'genes_mt' in file:
                f = open(top_dirs[id] + '/' + file)
                for i, line in enumerate(f):
                    if i == top_indices[id]:
                        top_genes_mt.append(line.replace('\n', '').split('\t'))
                f.close()
elif experiment_type == 'nuc':
    for id in range(0, len(top_dirs)):
        for file in os.listdir(top_dirs[id]):
            if 'genes_nuc' in file:
                f = open(top_dirs[id] + '/' + file)
                for i, line in enumerate(f):
                    if i == top_indices[id]:
                        top_genes_nuc.append(line.replace('\n', '').split('\t'))
                f.close()
else:
    for id in range(0, len(top_dirs)):
        for file in os.listdir(top_dirs[id]):
            if 'genes_nuc' in file:
                f = open(top_dirs[id] + '/' + file)
                for i, line in enumerate(f):
                    if i == top_indices[id]:
                        top_genes_nuc.append(line.replace('\n', '').split('\t'))
                f.close()
            if 'genes_mt' in file:
                f = open(top_dirs[id] + '/' + file)
                for i, line in enumerate(f):
                    if i == top_indices[id]:
                        top_genes_mt.append(line.replace('\n', '').split('\t'))
                f.close()

if len(result_file_suffix) > 0:
    result_file_suffix = '_' + result_file_suffix

if experiment_type == 'mt':
    for mt_genes_set in top_genes_mt:
        genes_mt_task = [mt_genes_set]
        genes_nuc_task = [[]]
        json_list = json.dumps([genes_mt_task, genes_nuc_task]).encode('utf-8')

        hash = hashlib.md5(json_list).hexdigest()

        root = result_path
        local_path = '/' + experiment_type + '/rf_type_' + str(
            random_forest_type) + '/ref_' + reference_pop + '_target_' + target_pop + '/top/' + hash + '/'
        fn_path = root + local_path
        pathlib.Path(fn_path).mkdir(parents=True, exist_ok=True)

        fn_test = str(target_accuracy) + '_accuracy' + result_file_suffix + '.txt'

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

            random_forest(fn_path)

elif experiment_type == 'nuc':
    for nuc_genes_set in top_genes_nuc:
        genes_mt_task = [[]]
        genes_nuc_task = [nuc_genes_set]
        json_list = json.dumps([genes_mt_task, genes_nuc_task]).encode('utf-8')

        hash = hashlib.md5(json_list).hexdigest()

        root = result_path
        local_path = '/' + experiment_type + '/rf_type_' + str(
            random_forest_type) + '/ref_' + reference_pop + '_target_' + target_pop + '/top/' + hash + '/'
        fn_path = root + local_path
        pathlib.Path(fn_path).mkdir(parents=True, exist_ok=True)

        fn_test = str(target_accuracy) + '_accuracy' + result_file_suffix + '.txt'

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

            random_forest(fn_path)

else:
    combinations = [[], []]
    for subset_id in range(0, len(top_genes_mt)):
        combinations[0].append(top_genes_mt[subset_id])
        combinations[1].append(top_genes_nuc[subset_id])
    for gene_set_id in range(0, len(combinations[0])):
        genes_mt_task = [combinations[0][gene_set_id]]
        genes_nuc_task = [combinations[1][gene_set_id]]
        json_list = json.dumps([genes_mt_task, genes_nuc_task]).encode('utf-8')

        hash = hashlib.md5(json_list).hexdigest()

        root = result_path
        local_path = '/' + experiment_type + '/rf_type_' + str(
            random_forest_type) + '/ref_' + reference_pop + '_target_' + target_pop + '/top/' + hash + '/'
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

            random_forest(fn_path)
