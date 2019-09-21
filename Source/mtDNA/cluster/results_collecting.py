import os
from Source.mtDNA.cluster.rf import random_forest

data_path = '../../../Data/'
data_path_npz = '../../../Data/genes/npz/'
data_path_pkl = '../../../Data/genes/pkl/'
experiment_type = 'nuc'
reference_pop = 'IBS'
target_pop = 'FIN'
reference_part = 0.75
result_file_suffix = 'cold'
target_accuracy = 0.68
num_features = 50
gene_files = ['test_gene_list_cold_adaptation.txt']
top_results_num = 5

result_path = '../../../Result/files/'
experiment_result_path = result_path + experiment_type + '/' + \
                         'ref(' + reference_pop + ')_' + 'target(' + target_pop + ')/'

top_accuracy = []
top_indices = []
top_dirs = []
top_genes_mt = []
top_genes_nuc = []

min_accuracy_value = 0.0
min_accuracy_index = top_results_num

dirnames = [x[0] for x in os.walk(experiment_result_path)][1:]

for dirname in dirnames:
    for filename in os.listdir(dirname):
        if 'accuracy' in filename:
            f = open(dirname + '/' + filename)
            line_num = 0
            for line in f:
                item = line.replace('\n', '')
                if len(top_accuracy) < top_results_num:
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
        k_mt_max = len(mt_genes_set)
        k_nuc_max = 1
        for k_mt in range(k_mt_max, k_mt_max + 1):
            for k_nuc in range(k_nuc_max, k_nuc_max + 1):
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
                    f.write('mt_genes_set\t' + ', '.join(mt_genes_set) + '\n')
                    f.write('k_mt\t' + str(k_mt) + '\n')
                    f.write('k_nuc\t' + str(k_nuc) + '\n')
                random_forest()

elif experiment_type == 'nuc':
    for nuc_genes_set in top_genes_nuc:
        k_mt_max = 1
        k_nuc_max = len(nuc_genes_set)
        for k_mt in range(k_mt_max, k_mt_max + 1):
            for k_nuc in range(k_nuc_max, k_nuc_max + 1):
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
                    f.write('nuc_genes_set\t' + ', '.join(nuc_genes_set) + '\n')
                    f.write('k_mt\t' + str(k_mt) + '\n')
                    f.write('k_nuc\t' + str(k_nuc) + '\n')
                random_forest()

else:
    for gene_set_id in range(0, len(top_genes_mt)):
        k_mt_max = len(top_genes_mt[gene_set_id])
        k_nuc_max = len(top_genes_nuc[gene_set_id])
        for k_mt in range(k_mt_max, k_mt_max + 1):
            for k_nuc in range(k_nuc_max, k_nuc_max + 1):
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
                    f.write('mt_genes_set\t' + ', '.join(top_genes_mt[gene_set_id]) + '\n')
                    f.write('nuc_genes_set\t' + ', '.join(top_genes_nuc[gene_set_id]) + '\n')
                    f.write('k_mt\t' + str(k_mt) + '\n')
                    f.write('k_nuc\t' + str(k_nuc) + '\n')
                random_forest()
