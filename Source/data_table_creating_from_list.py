import os

data_path = '../Data/'
data_chr_test_table_file_name = 'data_snp_test_genes_short.txt'
data_file = open(data_path + data_chr_test_table_file_name, 'w')

data_gene_list_file_name = 'test_gene_list_short.txt'
data_gene_file = open(data_path + data_gene_list_file_name)
genes = [line.replace('\n', '') for line in data_gene_file]
data_gene_file.close()

header = ''
for dir_name in os.listdir(data_path):
    if dir_name.startswith('chr'):
        if not dir_name.startswith('chrMT'):
            print(dir_name)
            chr_path = data_path + dir_name + '/'
            for gene_file_name in os.listdir(chr_path):
                f = open(chr_path + gene_file_name)
                for line in f:
                    if header == '':
                        header = line
                        data_file.write(header)
                    elif line == header:
                        continue
                    else:
                        curr_gene_data = line.split(' ')
                        gene_name = curr_gene_data[1]
                        if gene_name in genes:
                            data_file.write(line)
                f.close()
data_file.close()
