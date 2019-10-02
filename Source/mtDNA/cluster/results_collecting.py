from Source.mtDNA.cluster.rf import random_forest
import itertools
import json
import hashlib
import pathlib
import os

data_path = 'D:/Aaron/Bio/mtDNA/Data/'
data_path_npz = 'D:/Aaron/Bio/mtDNA/Data/genes/npz/'
data_path_pkl = 'D:/Aaron/Bio/mtDNA/Data/genes/pkl/'
experiment_type = 'mt-nuc'
reference_pop = 'IBS'
target_pop = 'FIN'
reference_part = 0.75
result_file_suffix = ''
target_accuracy = 0.6
num_features = 0
gene_files = ['mt_gene_list.txt', 'test_gene_list_short.txt']
create_tree = 0
k_mt_max = 1
k_nuc_max = 1
num_cluster_tasks = 6
num_atomic_tasks = 10
num_running_tasks = 0
num_top_results = 5

result_path = '../../../Result/files/'
experiment_result_path = result_path + experiment_type + '/' + \
                         'ref(' + reference_pop + ')_' + 'target(' + target_pop + ')/'

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

if experiment_type == 'mt':
    for mt_genes_set in top_genes_mt:
        genes_mt_task = mt_genes_set
        genes_nuc_task = []
        json_list = json.dumps([genes_mt_task, genes_nuc_task]).encode('utf-8')

        hash = hashlib.md5(json_list).hexdigest()

        root = 'D:/Aaron/Bio/mtDNA/Result/files'
        local_path = '/' + experiment_type + '/ref_' + reference_pop + '_target_' + target_pop + '/' + hash + '/'
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
        file_config.close()

        random_forest(fn_path)

elif experiment_type == 'nuc':
    for nuc_genes_set in top_genes_nuc:
        genes_mt_task = []
        genes_nuc_task = nuc_genes_set
        json_list = json.dumps([genes_mt_task, genes_nuc_task]).encode('utf-8')

        hash = hashlib.md5(json_list).hexdigest()

        root = 'D:/Aaron/Bio/mtDNA/Result/files'
        local_path = '/' + experiment_type + '/ref_' + reference_pop + '_target_' + target_pop + '/' + hash + '/'
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
        file_config.close()

        random_forest(fn_path)

else:
    combinations = [[], []]
    for subset_mt in top_genes_mt:
        for subset_nuc in top_genes_nuc:
            combinations[0].append(list(subset_mt))
            combinations[1].append(list(subset_nuc))
    for gene_set in combinations:
        genes_mt_task = gene_set[0]
        genes_nuc_task = gene_set[1]
        json_list = json.dumps([genes_mt_task, genes_nuc_task]).encode('utf-8')

        hash = hashlib.md5(json_list).hexdigest()

        root = 'D:/Aaron/Bio/mtDNA/Result/files'
        local_path = '/' + experiment_type + '/ref_' + reference_pop + '_target_' + target_pop + '/' + hash + '/'
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
        file_config.close()

        random_forest(fn_path)
