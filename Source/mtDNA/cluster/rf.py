from mtDNA.cluster.rf_unit import unit_task
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

    gene_id = 0
    for gene_list in genes:
        for item_id in range(0, len(gene_list)):
            item = gene_list[item_id]
            data_npz = np.load(data_path_npz + item + '.npz')
            unit_config.data.append(data_npz['data'])
            unit_config.data_position_dict[item] = gene_id
            gene_id += 1

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
            unit_config.person_index_dict.append(pickle.load(handle))
        with open(data_path_pkl + 'person_index_nuc_dict.pickle', 'rb') as handle:
            unit_config.person_index_dict.append(pickle.load(handle))

    with open(data_path_pkl + 'pop_person_dict.pickle', 'rb') as handle:
        unit_config.pop_person_dict = pickle.load(handle)

    results = Result()
    unit_task(unit_config, results)

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

    suffix = config_dict['result_file_suffix']
    if len(suffix) > 0:
        suffix += '_'

    if len(results.accuracy) > 0:
        with open(experiment_result_path + suffix + 'accuracy.txt', 'w') as f:
            for item in results.accuracy:
                f.write("%s\n" % item)

        with open(experiment_result_path + suffix + 'num_features.txt', 'w') as f:
            for item in results.num_features:
                f.write("%s\n" % item)

        if unit_config.params_dict['experiment_type'] == 'mt':
            with open(experiment_result_path + suffix + 'genes_mt.txt', 'w') as f:
                for item in results.mt_genes:
                    f.write("%s\n" % '\t'.join(list(map(str, item))))
        elif unit_config.params_dict['experiment_type'] == 'nuc':
            with open(experiment_result_path + suffix + 'genes_nuc.txt', 'w') as f:
                for item in results.nuc_genes:
                    f.write("%s\n" % '\t'.join(list(map(str, item))))
        else:
            with open(experiment_result_path + suffix + 'genes_mt.txt', 'w') as f:
                for item in results.mt_genes:
                    f.write("%s\n" % '\t'.join(list(map(str, item))))
            with open(experiment_result_path + suffix + 'genes_nuc.txt', 'w') as f:
                for item in results.nuc_genes:
                    f.write("%s\n" % '\t'.join(list(map(str, item))))
