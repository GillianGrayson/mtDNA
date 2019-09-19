import numpy as np
import pickle

data_path = '../../Data/genes/'
data_path_npz = '../../Data/genes/npz/'
gene = 'MT-ND1'

data_npz = np.load(data_path_npz + gene + '.npz')
data = data_npz['data']
with open(data_path + 'gene_snp_dict.pickle', 'rb') as handle:
    data_pkl = pickle.load(handle)
with open(data_path + 'person_index_nuc_dict.pickle', 'rb') as handle:
    data_pkl = pickle.load(handle)
with open(data_path + 'person_index_mt_dict.pickle', 'rb') as handle:
    data_pkl = pickle.load(handle)
with open(data_path + 'gene_chr_dict.pickle', 'rb') as handle:
    data_pkl = pickle.load(handle)
olo = 9