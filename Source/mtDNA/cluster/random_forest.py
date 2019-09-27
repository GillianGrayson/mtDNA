from rf_unit import unit_task
import numpy as np
import pickle
import os, errno

class Config:
    def __init__(self):
        self.data = []
        self.params_dict = {}
        self.gene_snp_dict = {}
        self.person_index_dict = {}
        self.data_position_dict = {}
        self.pop_person_dict = {}

class Result:
    def __init__(self):
        self.accuracy = []
        self.mt_genes = []
        self.nuc_genes = []
        self.num_features = []
        self.features = {}


def random_forest(config_dict):

    result_path = ''

    genes = []
    genes_ids = []
    genes_set_feature = []

    if 'mt_genes_set_feature' in config_dict:
        genes_set_feature.append([int(i) for i in config_dict['mt_genes_set_feature']])
    if 'nuc_genes_set_feature' in config_dict:
        genes_set_feature.append([int(i) for i in config_dict['nuc_genes_set_feature']])

    for file_id in range(0, len(config_dict['gene_files'])):
        data_gene_file = open(config_dict['data_path'] + config_dict['gene_files'][file_id])
        if len(genes_set_feature) > 0:
            genes.append([line.replace('\n', '') for i, line in enumerate(data_gene_file) if
                          i in genes_set_feature[file_id]])
        else:
            genes.append([line.replace('\n', '') for i, line in enumerate(data_gene_file)])
        genes_ids.append([])
        for gene_id in range(0, len(genes[file_id])):
            genes_ids[file_id].append(gene_id)
        data_gene_file.close()

    config_dict['genes_list'] = genes
    config_dict['genes_ids_list'] = genes_ids

    unit_config = Config()

    gene_id = 0
    for gene_list in genes:
        for item_id in range(0, len(gene_list)):
            item = gene_list[item_id]
            data_npz = np.load(config_dict['data_path_npz'] + item + '.npz')
            unit_config.data.append(data_npz['data'])
            unit_config.data_position_dict[item] = gene_id
            gene_id += 1

    unit_config.params_dict = config_dict

    with open(config_dict['data_path_pkl'] + 'gene_snp_dict.pickle', 'rb') as handle:
        unit_config.gene_snp_dict = pickle.load(handle)

    if unit_config.params_dict['experiment_type'] == 'mt':
        with open(config_dict['data_path_pkl'] + 'person_index_mt_dict.pickle', 'rb') as handle:
            unit_config.person_index_dict = pickle.load(handle)
    elif unit_config.params_dict['experiment_type'] == 'nuc':
        with open(config_dict['data_path_pkl'] + 'person_index_nuc_dict.pickle', 'rb') as handle:
            unit_config.person_index_dict = pickle.load(handle)
    else:
        unit_config.person_index_dict = []
        with open(config_dict['data_path_pkl'] + 'person_index_mt_dict.pickle', 'rb') as handle:
            unit_config.person_index_dict.append(pickle.load(handle))
        with open(config_dict['data_path_pkl'] + 'person_index_nuc_dict.pickle', 'rb') as handle:
            unit_config.person_index_dict.append(pickle.load(handle))

    with open(config_dict['data_path_pkl'] + 'pop_person_dict.pickle', 'rb') as handle:
        unit_config.pop_person_dict = pickle.load(handle)

    results = Result()
    unit_task(unit_config, results)

    if int(config_dict['create_tree']) == 1:
        experiment_result_path = result_path + \
                                 unit_config.params_dict['experiment_type'] + '/' + \
                                 'ref(' + unit_config.params_dict['reference_pop'] + ')_' + \
                                 'target(' + unit_config.params_dict['target_pop'] + ')/' + \
                                 'mt(' + str(unit_config.params_dict['k_mt']) + ')_' + \
                                 'nuc(' + str(unit_config.params_dict['k_nuc']) + ')/'

        try:
            os.makedirs(experiment_result_path)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
    else:
        experiment_result_path = result_path

    suffix = config_dict['result_file_suffix']
    if len(suffix) > 0:
        suffix = '_' + suffix

    if int(unit_config.params_dict['num_features']) == 0:
        with open(experiment_result_path + str(config_dict['target_accuracy']) +
                  '_accuracy' + suffix + '.txt', 'w') as f:
            for item in results.accuracy:
                f.write("%s\n" % item)

        with open(experiment_result_path + str(config_dict['target_accuracy']) +
                  '_num_features' + suffix + '.txt', 'w') as f:
            for item in results.num_features:
                f.write("%s\n" % item)

        if unit_config.params_dict['experiment_type'] == 'mt':
            with open(experiment_result_path + str(config_dict['target_accuracy']) +
                      '_genes_mt' + suffix + '.txt', 'w') as f:
                for item in results.mt_genes:
                    f.write("%s\n" % '\t'.join(list(map(str, item))))
        elif unit_config.params_dict['experiment_type'] == 'nuc':
            with open(experiment_result_path + str(config_dict['target_accuracy']) +
                      '_genes_nuc' + suffix + '.txt', 'w') as f:
                for item in results.nuc_genes:
                    f.write("%s\n" % '\t'.join(list(map(str, item))))
        else:
            with open(experiment_result_path + str(config_dict['target_accuracy']) +
                      '_genes_mt' + suffix + '.txt', 'w') as f:
                for item in results.mt_genes:
                    f.write("%s\n" % '\t'.join(list(map(str, item))))
            with open(experiment_result_path + str(config_dict['target_accuracy']) +
                      '_genes_nuc' + suffix + '.txt', 'w') as f:
                for item in results.nuc_genes:
                    f.write("%s\n" % '\t'.join(list(map(str, item))))
    else:
        if unit_config.params_dict['experiment_type'] == 'mt':
            with open(experiment_result_path + str(config_dict['target_accuracy']) +
                      '_top_features_mt(' + '_'.join(config_dict['mt_genes_set_feature']) + ')' +
                      suffix + '.txt', 'w') as f:
                features_count = 0
                for key, value in results.features.items():
                    if value > 0.0 and features_count < int(unit_config.params_dict['num_features']):
                        line = str(key) + '\t' + str(value)
                        f.write(line)
                        f.write('\n')
                        features_count += 1
        elif unit_config.params_dict['experiment_type'] == 'nuc':
            with open(experiment_result_path + str(config_dict['target_accuracy']) +
                      '_top_features_nuc(' + '_'.join(config_dict['nuc_genes_set_feature']) + ')' +
                      suffix + '.txt', 'w') as f:
                features_count = 0
                for key, value in results.features.items():
                    if value > 0.0 and features_count < int(unit_config.params_dict['num_features']):
                        line = str(key) + '\t' + str(value)
                        f.write(line)
                        f.write('\n')
                        features_count += 1
        else:
            with open(experiment_result_path + str(config_dict['target_accuracy']) +
                      '_top_features_mt(' + '_'.join(config_dict['mt_genes_set_feature']) +
                      ')_nuc(' + '_'.join(config_dict['nuc_genes_set_feature']) + ')' +
                      suffix + '.txt', 'w') as f:
                features_count = 0
                for key, value in results.features.items():
                    if value > 0.0 and features_count < int(unit_config.params_dict['num_features']):
                        line = str(key) + '\t' + str(value)
                        f.write(line)
                        f.write('\n')
                        features_count += 1


config_dict = {}
f = open('config.txt')
for line in f:
    line = line.replace('\n', '')
    items = line.split('\t')
    if items[0] == 'gene_files' or items[0] == 'mt_genes_set_feature' or items[0] == 'nuc_genes_set_feature':
        config_dict[items[0]] = items[1].split(', ')
    else:
        config_dict[items[0]] = items[1]
f.close()

config_dict['config_mt_genes'] = []
config_dict['config_nuc_genes'] = []
config_dict['config_mt_genes_set'] = []
config_dict['config_nuc_genes_set'] = []

f = open('config_mt_genes.txt')
for line in f:
    line = line.replace('\n', '')
    if len(line) > 0:
        items = line.split('\t')
        config_dict['config_mt_genes'].append([int(item) for item in items])
f.close()

f = open('config_nuc_genes.txt')
for line in f:
    line = line.replace('\n', '')
    if len(line) > 0:
        items = line.split('\t')
        config_dict['config_nuc_genes'].append([int(item) for item in items])
f.close()
config_dict['config_mt_genes_set'] = list(
    set([item for sublist in config_dict['config_mt_genes'] for item in sublist]))
config_dict['config_nuc_genes_set'] = list(
    set([item for sublist in config_dict['config_nuc_genes'] for item in sublist]))
random_forest(config_dict)