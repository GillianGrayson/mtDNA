import os

data_path = '../Data/'
data_table_file_name = 'data_snp.txt'
data_file = open(data_path + data_table_file_name, 'w')

header = ''
for dir_name in os.listdir(data_path):
    if dir_name.startswith('chr'):
        print(dir_name)
        chr_path = data_path + dir_name + '/'
        for gene_file_name in os.listdir(chr_path):
            f = open(chr_path + gene_file_name)
            for line in f:
                if header == '':
                    header = line
                    data_file.write(header)
                else:
                    data_file.write(line)
            f.close()
data_file.close()
