import os

data_path = '../../../Data/'
experiment_type = 'nuc'
reference_pop = 'IBS'
target_pop = 'FIN'
reference_part = 0.75
result_file_suffix = 'cold'
target_accuracy = 0.68
num_features = 0
gene_files = ['test_gene_list_cold_adaptation.txt']
k_mt_max = 1
k_nuc_max = 2
top_results_num = 5

result_path = '../../../Result/files/'
experiment_result_path = result_path + experiment_type + '/' + \
                         'ref(' + reference_pop + ')_' + 'target(' + target_pop + ')/'

top_accuracy = []
top_genes = []

min_accuracy_value = 0.0
min_accuracy_index = top_results_num

for (dirpath, dirnames, filenames) in os.walk(experiment_result_path):
    for filename in filenames:
        if filename.startswith('accuracy'):
            f = open(filename)
            for line in f:
                item = line.replace('\n', '')
                if float(item) > min_accuracy_value: