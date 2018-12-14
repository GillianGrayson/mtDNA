data_path = '../Data/'
data_file_name = 'data_snp_test.txt'

pop_file_name = 'samples_populations.txt'
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

target_chromosomes = ['chr22','chrMT']
snp_nuclear_variations = ['0|0', '0|1', '1|0', '1|1']
snp_mt_variations = ['0', '1']

snps_chrs = []
snps_names = []
snps_pos = []
snps_freqs = []
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
        snp_freqs = {}
        if snp_chr == 'chrMT':
            for id in range(0, len(snp_data)):
                sample_name = samples_names[id]
                for k in pop_dict.keys():
                    if sample_name in pop_dict[k]:
                        snp_pop = k
                if snp_pop in list(snp_freqs.keys()):
                    snp_freqs[snp_pop][snp_data[id]] += 1
                else:
                    snp_freqs[snp_pop] = dict.fromkeys(snp_mt_variations, 0)
                    snp_freqs[snp_pop][snp_data[id]] += 1
        else:
            for id in range(0, len(snp_data)):
                sample_name = samples_names[id]
                for k in pop_dict.keys():
                    if sample_name in pop_dict[k]:
                        snp_pop = k
                if snp_pop in list(snp_freqs.keys()):
                    snp_freqs[snp_pop][snp_data[id]] += 1
                else:
                    snp_freqs[snp_pop] = dict.fromkeys(snp_nuclear_variations, 0)
                    snp_freqs[snp_pop][snp_data[id]] += 1
        snps_chrs.append(snp_chr)
        snps_names.append(snp_name)
        snps_pos.append(snp_pos)
        snps_freqs.append(snp_freqs)

f.close()


