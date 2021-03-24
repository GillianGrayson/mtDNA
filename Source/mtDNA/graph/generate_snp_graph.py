import os
import pickle
import numpy as np
from tqdm import tqdm

root_path = 'E:/YandexDisk/mtDNA/'
data_path = root_path + 'Data/'
data_path_npz = data_path + 'genes/npz/'
data_path_pkl = data_path + 'genes/pkl/'

result_path = root_path + 'Result/graphs/snp/'
if not os.path.exists(result_path):
    os.makedirs(result_path)

result_path_files = result_path + 'files/'
if not os.path.exists(result_path_files):
    os.makedirs(result_path_files)

gene_list = []
f = open(data_path + 'test_gene_list_cold_adaptation.txt')
for line in f:
    gene_list.append(line.rstrip())
f.close()
f = open(data_path + 'mt_gene_list.txt')
for line in f:
    gene_list.append(line.rstrip())
f.close()

with open(data_path_pkl + 'person_index_mt_dict.pickle', 'rb') as handle:
    person_index_mt_dict = pickle.load(handle)

with open(data_path_pkl + 'person_index_nuc_dict.pickle', 'rb') as handle:
    person_index_nuc_dict = pickle.load(handle)

with open(data_path_pkl + 'gene_snp_dict.pickle', 'rb') as handle:
    gene_snp_dict = pickle.load(handle)

common_subjects = []
for subject in person_index_mt_dict:
    if subject in person_index_nuc_dict:
        common_subjects.append(subject)

nuc_mutation_dict = {subject: {} for subject in common_subjects}
mt_mutation_dict = {subject: {} for subject in common_subjects}
for gene_id in tqdm(range(len(gene_list))):
    gene = gene_list[gene_id]
    data_npz = np.load(data_path_npz + gene + '.npz')
    curr_gene_data = data_npz['data']
    snp_names = gene_snp_dict[gene]
    for i in range(0, len(common_subjects)):
        curr_subject = common_subjects[i]
        if gene.startswith('MT'):
            person_index = person_index_mt_dict[curr_subject]
            curr_subject_gene_data = curr_gene_data[:, person_index]
            for snp_name in snp_names:
                snp_id = snp_names[snp_name]
                if curr_subject_gene_data[snp_id] > 0:
                    mt_mutation_dict[curr_subject][snp_name] = 1
                else:
                    mt_mutation_dict[curr_subject][snp_name] = 0
        else:
            person_index = person_index_nuc_dict[curr_subject]
            curr_subject_gene_data = curr_gene_data[:, person_index]
            for snp_name in snp_names:
                snp_id = snp_names[snp_name]
                if curr_subject_gene_data[snp_id] > 0:
                    nuc_mutation_dict[curr_subject][snp_name] = 1
                else:
                    nuc_mutation_dict[curr_subject][snp_name] = 0

for subject in common_subjects:
    f = open(result_path_files + subject + '.txt', 'w+')
    mt_mutated_snps = [snp for snp in mt_mutation_dict[subject] if mt_mutation_dict[subject][snp] == 1]
    for nuc_snp in nuc_mutation_dict[subject]:
        if nuc_mutation_dict[subject][nuc_snp] == 1:
            f.write(nuc_snp + ':' + ','.join(mt_mutated_snps) + '\n')
        else:
            f.write(nuc_snp + ':' + '\n')
    nuc_mutated_snps = [snp for snp in nuc_mutation_dict[subject] if nuc_mutation_dict[subject][snp] == 1]
    for mt_snp in mt_mutation_dict[subject]:
        if mt_mutation_dict[subject][mt_snp] == 1:
            f.write(mt_snp + ':' + ','.join(nuc_mutated_snps) + '\n')
        else:
            f.write(mt_snp + ':' + '\n')
    f.close()

f = open(result_path + 'gene_snp.txt', 'w+')
f.write('SNP\tGene\n')
for gene in gene_snp_dict:
    if gene in gene_list:
        for snp in gene_snp_dict[gene]:
            f.write(snp + '\t' + gene + '\n')
f.close()
