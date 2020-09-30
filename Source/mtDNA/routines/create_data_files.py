import os
import numpy as np
import pickle

data_path = 'E:/YandexDisk/mtDNA/Data/genes/pkl/'
data_path_txt = 'E:/YandexDisk/mtDNA/Data/chr/'
data_path_npz = 'E:/YandexDisk/mtDNA/Data/genes/npz/'

gene_snp_dict = dict()
person_index_nuc_dict = dict()
person_index_mt_dict = dict()
gene_chr_dict = dict()
pop_person_dict = dict()

exclude_chromosomes = ['X', 'Y']

for folder_name in os.listdir(data_path_txt):
    print('Chr ' + folder_name)
    if folder_name not in exclude_chromosomes:
        for file_name in os.listdir(data_path_txt + folder_name):
            print('Chr ' + folder_name + ': ' + file_name[:-4])
            with open(data_path_txt + folder_name + '/' + file_name) as f:
                for i, l in enumerate(f):
                    pass
            num_lines = i
            f = open(data_path_txt + folder_name + '/' + file_name)
            header = f.readline()
            header = header.replace('\n', '')
            subjects = header.split(' ')[15:]
            if len(subjects) == 0:
                subjects = header.split('\t')[15:]
            l = l.replace('\n', '')
            curr_data = l.split(' ')
            if len(curr_data) == 1:
                curr_data = l.split('\t')
            chr = curr_data[0]
            gene = curr_data[1]
            for i in range(0, len(subjects)):
                person = subjects[i]
                if chr == 'chrMT' or chr == 'MT':
                    if person in person_index_mt_dict:
                        if person_index_mt_dict[person] == i:
                            continue
                        else:
                            print('For gene ', gene, ', line ', str(i), 'person mismatch')
                    else:
                        person_index_mt_dict[person] = i
                else:
                    if person in person_index_nuc_dict:
                        if person_index_nuc_dict[person] == i:
                            continue
                        else:
                            print('For gene ', gene, ', line ', str(i), 'person mismatch')
                    else:
                        person_index_nuc_dict[person] = i
            gene_snp_dict[gene] = dict()
            data = np.empty(shape=(num_lines, len(subjects)), dtype=int)
            line_id = 0
            for line in f:
                line = line.replace('\n', '')
                curr_data = line.split(' ')
                if len(curr_data) == 1:
                    curr_data = l.split('\t')
                chr = curr_data[0]
                if chr.startswith('chr'):
                    chr = chr[3:]
                gene = curr_data[1]
                pos = curr_data[2]
                snp = curr_data[3]
                if chr == 'MT':
                    gene_snp_dict[gene][pos] = line_id
                else:
                    gene_snp_dict[gene][snp] = line_id
                if gene not in gene_chr_dict:
                    gene_chr_dict[gene] = chr
                subject_id = 0
                for item in curr_data[15:]:
                    if chr == 'MT':
                        if item == '0':
                            data[line_id, subject_id] = 0
                        elif item == '1':
                            data[line_id, subject_id] = 1
                    else:
                        if item == '0|0':
                            data[line_id, subject_id] = 0
                        elif item == '0|1':
                            data[line_id, subject_id] = 1
                        elif item == '1|0':
                            data[line_id, subject_id] = 2
                        elif item == '1|1':
                            data[line_id, subject_id] = 3
                    subject_id += 1
                line_id += 1
            f.close()
            np.savez_compressed(data_path_npz + gene, data=data)

pop_file_name = 's_pop.txt'
f = open(data_path[:-10] + pop_file_name)
f.readline()
for line in f:
    line = line.replace('\n', '')
    curr_pop_data = line.split('\t')
    if curr_pop_data[1] in pop_person_dict:
        pop_person_dict[curr_pop_data[1]].append(curr_pop_data[0])
    else:
        pop_person_dict[curr_pop_data[1]] = []
        pop_person_dict[curr_pop_data[1]].append(curr_pop_data[0])
f.close()

with open(data_path + 'gene_snp_dict.pickle', 'wb') as handle:
    pickle.dump(gene_snp_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
with open(data_path + 'person_index_nuc_dict.pickle', 'wb') as handle:
    pickle.dump(person_index_nuc_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
with open(data_path + 'person_index_mt_dict.pickle', 'wb') as handle:
    pickle.dump(person_index_mt_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
with open(data_path + 'gene_chr_dict.pickle', 'wb') as handle:
    pickle.dump(gene_chr_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
with open(data_path + 'pop_person_dict.pickle', 'wb') as handle:
    pickle.dump(pop_person_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
