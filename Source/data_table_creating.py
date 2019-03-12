import os

data_path = '../Data/'
data_chr_table_file_name = 'data_snp.txt'
data_mt_table_file_name = 'data_snp_mt.txt'
data_file = open(data_path + data_chr_table_file_name, 'w')
data_mt_file = open(data_path + data_mt_table_file_name, 'w')

header = ''
for dir_name in os.listdir(data_path):
    if dir_name.startswith('chr'):
        if dir_name.startswith('chrMT'):
            header = ''
            chr_path = data_path + dir_name + '/'
            for gene_file_name in os.listdir(chr_path):
                f = open(chr_path + gene_file_name)
                for line in f:
                    if header == '':
                        header = line
                        data_mt_file.write(header)
                    elif line == header:
                        continue
                    else:
                        data_mt_file.write(line)
                f.close()
        else:
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
                        data_file.write(line)
                f.close()
data_file.close()
data_mt_file.close()
