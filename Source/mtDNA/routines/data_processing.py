import os

data_path = 'E:/YandexDisk/mtDNA/Data/'
genes_file_name = 'mart_export.txt'
target_chromosome = 'X'
result_path = 'E:/YandexDisk/mtDNA/Data/chr/' + target_chromosome + '/'
if not os.path.exists(result_path):
    os.makedirs(result_path)

f = open(data_path + genes_file_name)
f.readline()
gene_data = []
for line in f:
    curr_gene_data = line.rstrip().split('\t')
    gene_name = curr_gene_data[0]
    chromosome_name = curr_gene_data[4]
    gene_start = int(curr_gene_data[2])
    gene_finish = int(curr_gene_data[3])
    if chromosome_name == target_chromosome:
        gene_data.append([chromosome_name, gene_name, gene_start, gene_finish])
f.close()

if target_chromosome == 'MT':
    file_name = 'ALL.chr' + target_chromosome + '.phase3_callmom-v0_4.20130502.genotypes.vcf'
elif target_chromosome == 'X':
    file_name = 'ALL.chr' + target_chromosome + '.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf'
elif target_chromosome == 'Y':
    file_name = 'ALL.chr' + target_chromosome + '.phase3_integrated_v2a.20130502.genotypes.vcf'
else:
    file_name = 'ALL.chr' + target_chromosome + '.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf'
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

            if len(snp_data) > 0:

                fn = result_path + gene_name + '.txt'
                snp_data[0] = column_data
                fw = open(fn, 'w')
                for item in snp_data:
                    fw.write('\t'.join(item) + '\n')
                fw.close()
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
            column_data.append(column_names[column_id])
            curr_data.append(row[column_id])
        if len(set(curr_data[6:])) == 1:
            continue
    else:
        continue

    snp_data.append(curr_data)
f.close()
