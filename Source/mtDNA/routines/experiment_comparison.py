import os
import json
import hashlib
import pandas as pd

data_path = 'E:/YandexDisk/mtDNA/Result/files/'
gene_path = 'E:/YandexDisk/mtDNA/Data/'
random_forest_type = 3
target_accuracy = 0.5
gene_files = ['mt_gene_list.txt', 'test_gene_list_cold_adaptation.txt']
reference_pops = ['GBR', 'FIN', 'TSI']
target_pops = ['GBR', 'FIN', 'TSI']
experiment_types = ['mt', 'nuc', 'mt-nuc']

mt_genes = {}
mt_features = {}
nuc_genes = {}
nuc_features = {}
mt_nuc_genes = {}
mt_nuc_features = {}
accuracy = {}
for reference in reference_pops:
    mt_genes[reference] = {}
    mt_features[reference] = {}
    nuc_genes[reference] = {}
    nuc_features[reference] = {}
    mt_nuc_genes[reference] = {}
    mt_nuc_features[reference] = {}
    accuracy[reference] = {}
    for target in target_pops:
        if target != reference:
            mt_genes[reference][target] = []
            mt_features[reference][target] = [[], []]
            nuc_genes[reference][target] = []
            nuc_features[reference][target] = [[], []]
            mt_nuc_genes[reference][target] = []
            mt_nuc_features[reference][target] = [[], []]
            accuracy[reference][target] = 0.0

for reference_pop in reference_pops:
    for target_pop in target_pops:
        if reference_pop != target_pop:
            genes_mt = []
            genes_names_mt = []
            genes_nuc = []
            genes_names_nuc = []
            for file_id in range(0, len(gene_files)):
                data_gene_file = open(gene_path + gene_files[file_id])
                if file_id == 0:
                    for i, line in enumerate(data_gene_file):
                        genes_mt.append(i)
                        genes_names_mt.append(line.rstrip())
                    genes_nuc = []
                    genes_names_nuc = []
                else:
                    for i, line in enumerate(data_gene_file):
                        genes_nuc.append(i)
                        genes_names_nuc.append(line.rstrip())
                data_gene_file.close()

            for experiment_type in experiment_types:

                if experiment_type == 'mt':
                    data_hash = [genes_mt, []]
                if experiment_type == 'nuc':
                    data_hash = [[], genes_nuc]
                if experiment_type == 'mt-nuc':
                    data_hash = [genes_mt, genes_nuc]

                json_list = json.dumps(data_hash).encode('utf-8')
                curr_hash = hashlib.md5(json_list).hexdigest()

                path = data_path + experiment_type + '/rf_type_' + str(random_forest_type) + \
                       '/ref_' + reference_pop + '_target_' + target_pop + '/' + curr_hash + '/'

                top_filename = str(target_accuracy) + '_top.txt'
                f = open(path + top_filename)
                curr_accuracy = float(f.readline().rstrip())
                accuracy[reference_pop][target_pop] = curr_accuracy

                if experiment_type == 'mt' or experiment_type == 'mt-nuc':
                    curr_mt_genes = f.readline()
                if experiment_type == 'nuc' or experiment_type == 'mt-nuc':
                    curr_nuc_genes = f.readline()
                curr_features = f.readline().rstrip()
                curr_features_list = curr_features.split(';')[:-1]
                for feature in curr_features_list:
                    feature_name = feature.split(':')[0]
                    feature_value = float(feature.split(':')[1])
                    if experiment_type == 'mt':
                        mt_features[reference_pop][target_pop][0].append(feature_name)
                        mt_features[reference_pop][target_pop][1].append(feature_value)
                    if experiment_type == 'nuc':
                        nuc_features[reference_pop][target_pop][0].append(feature_name)
                        nuc_features[reference_pop][target_pop][1].append(feature_value)
                    if experiment_type == 'mt-nuc':
                        mt_nuc_features[reference_pop][target_pop][0].append(feature_name)
                        mt_nuc_features[reference_pop][target_pop][1].append(feature_value)
                f.close()

result_dict = {'reference': [], 'target': [], 'mt-nuc': [], 'mt-nuc_importance': [], 'mt-nuc_importance_norm': [],
               'mt_change': [], 'mt_change_norm': [], 'nuc_change': [], 'nuc_change_norm': []}
for reference_pop in reference_pops:
    for target_pop in target_pops:
        if reference_pop != target_pop:
            for feature_id in range(0, len(mt_nuc_features[reference_pop][target_pop][0])):
                feature_name = mt_nuc_features[reference_pop][target_pop][0][feature_id]
                mt_feature = feature_name.split('_')[0]
                nuc_feature = feature_name.split('_')[1]
                importance = mt_nuc_features[reference_pop][target_pop][1][feature_id]
                importance_min = min(mt_nuc_features[reference_pop][target_pop][1])
                importance_max = max(mt_nuc_features[reference_pop][target_pop][1])
                importance_norm = (importance - importance_min) / (importance_max - importance_min)

                if mt_feature in mt_features[reference_pop][target_pop][0]:
                    mt_id = mt_features[reference_pop][target_pop][0].index(mt_feature)
                    mt_importance = mt_features[reference_pop][target_pop][1][mt_id]
                    mt_change = importance - mt_importance
                    mt_importance_min = min(mt_features[reference_pop][target_pop][1])
                    mt_importance_max = max(mt_features[reference_pop][target_pop][1])
                    mt_importance_norm = (mt_importance - mt_importance_min) / (mt_importance_max - mt_importance_min)
                    mt_change_norm = importance_norm - mt_importance_norm
                else:
                    mt_change = ''
                    mt_change_norm = ''

                if nuc_feature in nuc_features[reference_pop][target_pop][0]:
                    nuc_id = nuc_features[reference_pop][target_pop][0].index(nuc_feature)
                    nuc_importance = nuc_features[reference_pop][target_pop][1][nuc_id]
                    nuc_change = importance_norm - nuc_importance
                    nuc_importance_min = min(nuc_features[reference_pop][target_pop][1])
                    nuc_importance_max = max(nuc_features[reference_pop][target_pop][1])
                    nuc_importance_norm = (nuc_importance - nuc_importance_min) / (
                                nuc_importance_max - nuc_importance_min)
                    nuc_change_norm = importance_norm - nuc_importance_norm
                else:
                    nuc_change = ''
                    nuc_change_norm = ''

                result_dict['reference'].append(reference_pop)
                result_dict['target'].append(target_pop)
                result_dict['mt-nuc'].append(feature_name)
                result_dict['mt-nuc_importance'].append(importance)
                result_dict['mt-nuc_importance_norm'].append(importance_norm)
                result_dict['mt_change'].append(mt_change)
                result_dict['mt_change_norm'].append(mt_change_norm)
                result_dict['nuc_change'].append(nuc_change)
                result_dict['nuc_change_norm'].append(nuc_change_norm)

save_path = data_path + 'comparison/'
if not os.path.exists(save_path):
    os.makedirs(save_path)

result_df = pd.DataFrame(result_dict)
writer = pd.ExcelWriter(save_path + 'GBR_FIN_TSI.xlsx', engine='xlsxwriter')
result_df.to_excel(writer, index=False)
writer.save()
