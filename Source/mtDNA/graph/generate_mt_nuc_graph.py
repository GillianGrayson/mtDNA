import os
import pickle
import numpy as np
import pandas as pd
from tqdm import tqdm

root_path = 'E:/YandexDisk/mtDNA/'
data_path = root_path + 'Data/'
data_path_npz = data_path + 'genes/npz/'
data_path_pkl = data_path + 'genes/pkl/'

result_path = root_path + 'Result/graphs/'
if not os.path.exists(result_path):
    os.makedirs(result_path)

result_path_files = result_path + 'files/'
if not os.path.exists(result_path_files):
    os.makedirs(result_path_files)

with open(data_path_pkl + 'gene_chr_dict.pickle', 'rb') as handle:
    gene_chr_dict = pickle.load(handle)

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
for gene_id in tqdm(range(len(list(gene_chr_dict.keys())))):
    gene = list(gene_chr_dict.keys())[gene_id]
    data_npz = np.load(data_path_npz + gene + '.npz')
    curr_gene_data = data_npz['data']
    for i in range(0, len(common_subjects)):
        curr_subject = common_subjects[i]
        if gene_chr_dict[gene] == 'MT':
            person_index = person_index_mt_dict[curr_subject]
            curr_subject_gene_data = curr_gene_data[:, person_index]
            if len(set(curr_subject_gene_data)) > 1:
                mt_gene_mutation_dict[curr_subject][gene] = 1
            elif len(set(curr_subject_gene_data)) == 1 and curr_subject_gene_data[0] == 1:
                mt_gene_mutation_dict[curr_subject][gene] = 1
            else:
                mt_gene_mutation_dict[curr_subject][gene] = 0
        else:
            person_index = person_index_nuc_dict[curr_subject]
            curr_subject_gene_data = curr_gene_data[:, person_index]
            if len(set(curr_subject_gene_data)) > 1:
                nuc_gene_mutation_dict[curr_subject][gene] = 1
            elif len(set(curr_subject_gene_data)) == 1 and curr_subject_gene_data[0] > 0:
                nuc_gene_mutation_dict[curr_subject][gene] = 1
            else:
                nuc_gene_mutation_dict[curr_subject][gene] = 0

with open(result_path + 'mt_gene_mutation.pickle', 'wb') as handle:
    pickle.dump(mt_gene_mutation_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
with open(result_path + 'nuc_gene_mutation.pickle', 'wb') as handle:
    pickle.dump(nuc_gene_mutation_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

for subject in common_subjects:
    f = open(result_path_files + subject + '.txt', 'w+')
    mt_mutated_genes = [gene for gene in mt_gene_mutation_dict[subject] if mt_gene_mutation_dict[subject][gene] == 1]
    for nuc_gene in nuc_gene_mutation_dict[subject]:
        if nuc_gene_mutation_dict[subject][nuc_gene] == 1:
            f.write(nuc_gene + ':' + ','.join(mt_mutated_genes) + '\n')
        else:
            f.write(nuc_gene + ':' + '\n')
    nuc_mutated_genes = [gene for gene in nuc_gene_mutation_dict[subject] if nuc_gene_mutation_dict[subject][gene] == 1]
    for mt_gene in mt_gene_mutation_dict[subject]:
        if mt_gene_mutation_dict[subject][mt_gene] == 1:
            f.write(mt_gene + ':' + ','.join(nuc_mutated_genes) + '\n')
        else:
            f.write(mt_gene + ':' + '\n')
    f.close()

population_data = pd.read_csv(data_path + 's_pop.txt', delimiter='\t').to_dict('list')
subject_info = {'subject': [], 'population': [], 'super_population': []}
for subject in common_subjects:
    subject_index = population_data['sample'].index(subject)
    pop = population_data['pop'][subject_index]
    super_pop = population_data['super_pop'][subject_index]
    subject_info['subject'].append(subject)
    subject_info['population'].append(pop)
    subject_info['super_population'].append(super_pop)
df_info = pd.DataFrame(subject_info)
df_info.to_csv(result_path + 'subjects_info.txt', header=True, index=False, sep='\t')
