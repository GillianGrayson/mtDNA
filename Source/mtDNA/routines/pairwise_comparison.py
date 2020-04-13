import pandas as pd

data_path = '../Data/'
data_file_name = 'data_snp_test_genes.txt'
data_file_name_mt = 'data_snp_mt.txt'

pop_file_name = 's_pop.txt'
pop_dict = {}
f = open(data_path + pop_file_name)
for line in f:
    line = line.replace('\n', '')
    curr_pop_data = line.split(' ')
    if curr_pop_data[1] in pop_dict:
        pop_dict[curr_pop_data[1]].append(curr_pop_data[0])
    else:
        pop_dict[curr_pop_data[1]] = [curr_pop_data[0]]
f.close()

snp_nuclear_variations = ['0|0', '0|1', '1|0', '1|1']
snp_mt_variations = ['0', '1']

snps_chrs_mt = []
snps_names_mt = []
snps_pos_mt = []
snps_stats_mt = []
f = open(data_path + data_file_name_mt)
header = f.readline().replace('\n', '')
header = header.split(' ')
samples_names = header[6:]
for line in f:
    line = line.replace('\n', '')
    curr_snp_data = line.split(' ')
    snp_chr = curr_snp_data[0]
    snp_gene = curr_snp_data[1]
    snp_pos = curr_snp_data[2]
    snp_name = curr_snp_data[3]
    snp_data = curr_snp_data[6:]

    if snp_chr == 'chrMT':
        snp_stats_mt = {}
        for id in range(0, len(snp_data)):
            sample_name = samples_names[id]
            for k in pop_dict.keys():
                if sample_name in pop_dict[k]:
                    snp_pop = k
            if snp_pop in list(snp_stats_mt.keys()):
                snp_stats_mt[snp_pop][snp_data[id]] += 1
            else:
                snp_stats_mt[snp_pop] = dict.fromkeys(snp_mt_variations, 0)
                snp_stats_mt[snp_pop][snp_data[id]] += 1
        snps_chrs_mt.append(snp_chr)
        snps_names_mt.append(snp_name)
        snps_pos_mt.append(snp_pos)
        snps_stats_mt.append(snp_stats_mt)
f.close()
for snp_id in range(0, len(snps_names_mt)):
    for pop in list(snps_stats_mt[snp_id].keys()):
        sum_snps = sum(list(snps_stats_mt[snp_id][pop].values()))
        for snp_type in list(snps_stats_mt[snp_id][pop].keys()):
            snps_stats_mt[snp_id][pop][snp_type] = float(snps_stats_mt[snp_id][pop][snp_type]) / float(sum_snps)

snps_chrs = []
snps_names = []
snps_pos = []
snps_stats = []
f = open(data_path + data_file_name)
header = f.readline().replace('\n', '')
header = header.split(' ')
samples_names = header[6:]
for line in f:
    line = line.replace('\n', '')
    curr_snp_data = line.split(' ')
    snp_chr = curr_snp_data[0]
    snp_gene = curr_snp_data[1]
    snp_pos = curr_snp_data[2]
    snp_name = curr_snp_data[3]
    snp_data = curr_snp_data[6:]

    if snp_chr in target_chromosomes:
        snp_stats = {}
        for id in range(0, len(snp_data)):
            sample_name = samples_names[id]
            for k in pop_dict.keys():
                if sample_name in pop_dict[k]:
                    snp_pop = k
            if snp_pop in list(snp_stats.keys()):
                snp_stats[snp_pop][snp_data[id]] += 1
            else:
                snp_stats[snp_pop] = dict.fromkeys(snp_nuclear_variations, 0)
                snp_stats[snp_pop][snp_data[id]] += 1
        snps_chrs.append(snp_chr)
        snps_names.append(snp_name)
        snps_pos.append(snp_pos)
        snps_stats.append(snp_stats)

f.close()

for snp_id in range(0, len(snps_names)):
    for pop in list(snps_stats[snp_id].keys()):
        sum_snps = sum(list(snps_stats[snp_id][pop].values()))
        for snp_type in list(snps_stats[snp_id][pop].keys()):
            snps_stats[snp_id][pop][snp_type] = float(snps_stats[snp_id][pop][snp_type]) / float(sum_snps)
