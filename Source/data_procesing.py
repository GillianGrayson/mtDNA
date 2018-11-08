import numpy as np
import os

data_path = '../Data/'
genes_file_name = 'genes_locations_sort.txt'
population_file_name = 's_pop.txt'
target_chromosome = 'chr7'
result_path = '../Data/' + target_chromosome + '/'
if not os.path.exists(result_path):
    os.makedirs(result_path)

f = open(data_path + genes_file_name)
gene_data = []
for line in f:
    curr_gene_data = line.split('\t')
    gene_name = curr_gene_data[0]
    chromosome_name = curr_gene_data[1]
    gene_start = int(curr_gene_data[2])
    gene_finish = int(curr_gene_data[3])
    if chromosome_name == target_chromosome:
        gene_data.append([chromosome_name, gene_name, gene_start, gene_finish])
f.close()

pop_dict = dict()
pop_dict['GBR'] = 'N_EUR'
pop_dict['FIN'] = 'N_EUR'
pop_dict['IBS'] = 'S_EUR'
pop_dict['TSI'] = 'S_EUR'
pop_dict['YRI'] = 'AFR'
pop_dict['LWK'] = 'AFR'
pop_dict['GWD'] = 'AFR'
pop_dict['MSL'] = 'AFR'
pop_dict['ESN'] = 'AFR'
f = open(data_path + population_file_name)
population_data = dict((x,[]) for x in set(pop_dict.values()))
for line in f:
    curr_population_data = line.split('\t')
    sample_name = curr_population_data[0]
    population_name = curr_population_data[1]
    if population_name in pop_dict.keys():
        pop_name = pop_dict[population_name]
        population_data[pop_name].append(sample_name)
f.close()

file_name = 'ALL.' + target_chromosome + '.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf'
f = open(data_path + file_name)
column_names = []
pop_names = []
snp_data = []
column_data = []
genes_passed = 0
for line in f:
    curr_data = []
    if line.startswith('##'):
        continue
    elif line.startswith('#CHROM'):
        line = line.replace('\n', '')
        column_names = line.split('\t')
        for column_id in range(0, len(column_names)):
            for key in list(population_data.keys()):
                if column_names[column_id] in population_data[key]:
                    pop_names.append([column_names[column_id], key])
        fn = data_path + 'samples_populations.txt'
        np.savetxt(fn, pop_names, fmt='%s')
    elif line.startswith(target_chromosome[3:]):
        line = line.replace('\n', '')
        row = line.split('\t')
        curr_position = int(row[1])
        curr_snp_name = row[2]
        curr_origin_snp = row[3]
        curr_mutate_snp = row[4]
        curr_snp_mutates = row[9:]

        # position check
        curr_gene_info = gene_data[genes_passed]
        chromosome_name = curr_gene_info[0]
        gene_name = curr_gene_info[1]
        gene_start = curr_gene_info[2]
        gene_end = curr_gene_info[3]
        if gene_start <= curr_position <= gene_end:
            curr_gene = gene_name
        elif curr_position > gene_end:

            fn = result_path + gene_name + '.txt'
            snp_data[0] = column_data
            np.savetxt(fn, snp_data, fmt='%s')
            snp_data = []

            genes_passed += 1
            print(genes_passed)

            if genes_passed < len(gene_data):
                curr_gene_info = gene_data[genes_passed]
                gene_name = curr_gene_info[1]
                gene_start = curr_gene_info[2]
                gene_end = curr_gene_info[3]
                if gene_start <= curr_position <= gene_end:
                    curr_gene = gene_name
                else:
                    continue
            else:
                break
        else:
            continue

        # snp length check
        if len(curr_origin_snp) > 1 or len(curr_mutate_snp) > 1:
            continue

        # variability check
        if len(set(curr_snp_mutates)) == 1:
            continue

        curr_data = [chromosome_name, gene_name, str(curr_position), curr_snp_name, curr_origin_snp, curr_mutate_snp]
        column_data = [column_names[0][1:], 'GENE', column_names[1], column_names[2], column_names[3], column_names[4]]

        for column_id in range(0, len(column_names)):
            for key in list(population_data.keys()):
                if column_names[column_id] in population_data[key]:
                    pop_names.append([column_names[column_id], key])
                    column_data.append(column_names[column_id])
                    curr_data.append(row[column_id])
        if len(set(curr_data[6:])) == 1:
            continue
    else:
        continue

    snp_data.append(curr_data)
f.close()
