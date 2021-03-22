import os
import pickle

data_path = 'E:/YandexDisk/mtDNA/Data/genes/pkl/'
data_path_txt = 'E:/YandexDisk/mtDNA/Data/genes/txt'

gene_snp_dict = dict()
person_index_nuc_dict = dict()
person_index_mt_dict = dict()

for file_name in os.listdir(data_path_txt):
    print(file_name[:-4])
    with open(data_path_txt + '/' + file_name) as f:
        for i, l in enumerate(f):
            pass
    num_lines = i
    f = open(data_path_txt + '/' + file_name)
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
    line_id = 0
    for line in f:
        line = line.replace('\n', '')
        curr_data = line.split(' ')
        if len(curr_data) == 1:
            curr_data = line.split('\t')
        chr = curr_data[0]
        if chr.startswith('chr'):
            chr = chr[3:]
        gene = curr_data[1]
        pos = curr_data[2]
        snp = curr_data[3]
        curr_snp_name = snp
        if chr == 'chrMT' or chr == 'MT':
            curr_snp_name = pos
        if len(set(curr_data[15:])) > 1:
            gene_snp_dict[gene][curr_snp_name] = line_id
            line_id += 1
    f.close()

if gene_snp_dict:
    with open(data_path + 'gene_snp_dict_cold.pickle', 'wb') as handle:
        pickle.dump(gene_snp_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
