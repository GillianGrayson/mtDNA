import os
import pickle
import numpy as np
from tqdm import tqdm

root_path = 'E:/YandexDisk/mtDNA/'
data_path = root_path + 'Data/'
data_path_txt = data_path + 'genes/txt/'
data_path_npz = data_path + 'genes/npz/'
data_path_pkl = data_path + 'genes/pkl/'

result_path = root_path + 'Result/graphs/weighted/'
if not os.path.exists(result_path):
    os.makedirs(result_path)

result_path_files = result_path + 'files/'
if not os.path.exists(result_path_files):
    os.makedirs(result_path_files)

gene_list = []
for filename in os.listdir(data_path_txt):
    gene_list.append(filename[:-4])

with open(data_path_pkl + 'person_index_mt_dict.pickle', 'rb') as handle:
    person_index_mt_dict = pickle.load(handle)

with open(data_path_pkl + 'person_index_nuc_dict.pickle', 'rb') as handle:
    person_index_nuc_dict = pickle.load(handle)

common_subjects = []
for subject in person_index_mt_dict:
    if subject in person_index_nuc_dict:
        common_subjects.append(subject)

nuc_gene_mutation_dict = {subject: {} for subject in common_subjects}
mt_gene_mutation_dict = {subject: {} for subject in common_subjects}
for gene_id in tqdm(range(len(gene_list))):
    gene = list(gene_list)[gene_id]
    data_npz = np.load(data_path_npz + gene + '.npz')
    curr_gene_data = data_npz['data']
    for i in range(0, len(common_subjects)):
        curr_subject = common_subjects[i]
        if gene.startswith('MT-'):
            person_index = person_index_mt_dict[curr_subject]
            curr_subject_gene_data = curr_gene_data[:, person_index]
            mt_gene_mutation_dict[curr_subject][gene] = np.count_nonzero(curr_subject_gene_data)
        else:
            person_index = person_index_nuc_dict[curr_subject]
            curr_subject_gene_data = curr_gene_data[:, person_index]
            nuc_gene_mutation_dict[curr_subject][gene] = np.count_nonzero(curr_subject_gene_data)

for subject in common_subjects:
    f = open(result_path_files + subject + '.txt', 'w+')
    mt_mutated_genes = [gene for gene in mt_gene_mutation_dict[subject] if mt_gene_mutation_dict[subject][gene] > 0]
    nuc_mutated_genes = [gene for gene in nuc_gene_mutation_dict[subject] if nuc_gene_mutation_dict[subject][gene] > 0]
    for nuc_gene in nuc_gene_mutation_dict[subject]:
        num_nuc_mutations = nuc_gene_mutation_dict[subject][nuc_gene]
        for mt_gene in mt_gene_mutation_dict[subject]:
            num_mt_mutations = mt_gene_mutation_dict[subject][mt_gene]
            f.write(nuc_gene + '\t' + mt_gene + '\t' + str(num_nuc_mutations * num_mt_mutations) + '\n')
    for mt_gene in mt_gene_mutation_dict[subject]:
        num_mt_mutations = mt_gene_mutation_dict[subject][mt_gene]
        for nuc_gene in nuc_gene_mutation_dict[subject]:
            num_nuc_mutations = nuc_gene_mutation_dict[subject][nuc_gene]
            f.write(mt_gene + '\t' + nuc_gene + '\t' + str(num_nuc_mutations * num_mt_mutations) + '\n')
    f.close()
