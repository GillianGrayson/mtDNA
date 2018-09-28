data_path = '../Data/'
genes_file_name = 'genes_locations_sort.txt'
chromosome = 'chr22'
f = open(data_path+genes_file_name)
gene_data = dict()
for line in f:
    curr_gene_data = line.split('\t')
    gene_name = curr_gene_data[0]
    chromosome_name = curr_gene_data[1]
    gene_start = int(curr_gene_data[2])
    gene_finish = int(curr_gene_data[3])
    if chromosome_name == chromosome:
        gene_data[gene_name] = [chromosome_name, gene_start, gene_finish]


file_name = 'ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf'
f = open(data_path+file_name)
column_names = []
for line in f:
    if line.startswith('##'):
        continue
    elif line.startswith('#CHROM'):
        column_names = line.split('\t')
    else:
        column = line.split('\t')
