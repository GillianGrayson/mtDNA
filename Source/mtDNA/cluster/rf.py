from mtDNA.cluster.rf_unit import unit_task
import numpy as np
import pickle

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
        self.accuracy_mt_genes = []
        self.accuracy_nuc_genes = []
        self.features = {}


def random_forest(k_mt, k_nuc):
    config_dict = {}
    f = open('config.txt')
    for line in f:
        line = line.replace('\n', '')
        items = line.split('\t')
        if items[0] == 'gene_files':
            config_dict[items[0]] = items[1].split(', ')
        else:
            config_dict[items[0]] = items[1]
    f.close()

    data_path = '../../../Data/'
    data_path_npz = '../../../Data/genes/npz/'
    data_path_pkl = '../../../Data/genes/pkl/'
    result_path = '../../../Result/files/'

    genes = []
    genes_ids = []

    for file_id in range(0, len(config_dict['gene_files'])):
        data_gene_file = open(data_path + config_dict['gene_files'][file_id])
        genes.append([line.replace('\n', '') for line in data_gene_file])
        genes_ids.append([])
        for gene_id in range(0, len(genes[file_id])):
            genes_ids[file_id].append(gene_id)
        data_gene_file.close()

    config_dict['genes_list'] = genes
    config_dict['genes_ids_list'] = genes_ids

    unit_config = Config()

    for gene_list in genes:
        for item_id in range(0, len(gene_list)):
            item = gene_list[item_id]
            data_npz = np.load(data_path_npz + item + '.npz')
            unit_config.data.append(data_npz['data'])
            unit_config.data_position_dict[item] = item_id

    unit_config.params_dict = config_dict
    unit_config.params_dict['k_mt'] = k_mt
    unit_config.params_dict['k_nuc'] = k_nuc

    with open(data_path_pkl + 'gene_snp_dict.pickle', 'rb') as handle:
        unit_config.gene_snp_dict = pickle.load(handle)

    if unit_config.params_dict['experiment_type'] == 'mt':
        with open(data_path_pkl + 'person_index_mt_dict.pickle', 'rb') as handle:
            unit_config.person_index_dict = pickle.load(handle)
    elif unit_config.params_dict['experiment_type'] == 'nuc':
        with open(data_path_pkl + 'person_index_nuc_dict.pickle', 'rb') as handle:
            unit_config.person_index_dict = pickle.load(handle)
    else:
        unit_config.person_index_dict = []
        with open(data_path_pkl + 'person_index_mt_dict.pickle', 'rb') as handle:
            unit_config.person_index_dict[0] = pickle.load(handle)
        with open(data_path_pkl + 'person_index_nuc_dict.pickle', 'rb') as handle:
            unit_config.person_index_dict[1] = pickle.load(handle)

    with open(data_path_pkl + 'pop_person_dict.pickle', 'rb') as handle:
        unit_config.pop_person_dict = pickle.load(handle)

    results = Result()
    unit_task(unit_config, results)

    with open(result_path + '/' +
              unit_config.params_dict['experiment_type'] + '/' +
              'ref(' + unit_config.params_dict['reference_pop'] + ')_' +
              'target(' + unit_config.params_dict['target_pop'] + ')/' +
              unit_config.params_dict['k_mt'] + '_' + unit_config.params_dict['k_nuc'] + '/' +
              'accuracy.txt', 'w') as f:
        for item in results.accuracy:
            f.write("%s\n" % item)

    if unit_config.params_dict['experiment_type'] == 'mt':
        with open(result_path + '/' +
                  unit_config.params_dict['experiment_type'] + '/' +
                  'ref(' + unit_config.params_dict['reference_pop'] + ')_' +
                  'target(' + unit_config.params_dict['target_pop'] + ')/' +
                  unit_config.params_dict['k_mt'] + '_' + unit_config.params_dict['k_nuc'] + '/' +
                  'accuracy_mt_genes.txt', 'w') as f:
            for item in results.accuracy_mt_genes:
                f.write('\t'.join(item) + '\n')
    elif unit_config.params_dict['experiment_type'] == 'nuc':
        with open(result_path + '/' +
                  unit_config.params_dict['experiment_type'] + '/' +
                  'ref(' + unit_config.params_dict['reference_pop'] + ')_' +
                  'target(' + unit_config.params_dict['target_pop'] + ')/' +
                  unit_config.params_dict['k_mt'] + '_' + unit_config.params_dict['k_nuc'] + '/' +
                  'accuracy_nuc_genes.txt', 'w') as f:
            for item in results.accuracy_nuc_genes:
                f.write('\t'.join(item) + '\n')
    else:
        with open(result_path + '/' +
                  unit_config.params_dict['experiment_type'] + '/' +
                  'ref(' + unit_config.params_dict['reference_pop'] + ')_' +
                  'target(' + unit_config.params_dict['target_pop'] + ')/' +
                  unit_config.params_dict['k_mt'] + '_' + unit_config.params_dict['k_nuc'] + '/' +
                  'accuracy_mt_genes.txt', 'w') as f:
            for item in results.accuracy_mt_genes:
                f.write('\t'.join(item) + '\n')
        with open(result_path + '/' +
                  unit_config.params_dict['experiment_type'] + '/' +
                  'ref(' + unit_config.params_dict['reference_pop'] + ')_' +
                  'target(' + unit_config.params_dict['target_pop'] + ')/' +
                  unit_config.params_dict['k_mt'] + '_' + unit_config.params_dict['k_nuc'] + '/' +
                  'accuracy_nuc_genes.txt', 'w') as f:
            for item in results.accuracy_nuc_genes:
                f.write('\t'.join(item) + '\n')
