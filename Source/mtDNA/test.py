import numpy as np


a = [0, 1, 2, 3, 4]
quartiles = np.percentile(a, [25, 75])
iqr = quartiles[1] - quartiles[0]
print(iqr)

# a = [0, 1, 2]
# b = np.asarray(a)
# b = b.reshape(3, 1)
# print(b)


# import pandas as pd
#
# sample_file = 'D:/Aaron/Bio/turin/test/GSE40279_sample_key.txt'
# sample_key = pd.read_csv(sample_file, sep="\t", dtype='str', header=None)
# print('Samples keys loaded')
#
# data = pd.read_csv('D:/Aaron/Bio/turin/test/GSE40279.txt', sep="\t", dtype='str')
# print('Data loaded')
# print('Size', data.shape)
# column_names = []
# column_names.append('ID_REF')
# for id in range(0, len(sample_key[2])):
#     if sample_key[1][id] in data.columns[id + 1]:
#         column_names.append(sample_key[2][id])
#     else:
#         print('error in ', id)
# data.columns = column_names
# data = data.set_index('ID_REF').T   #ID_REF for GSE40279, id for EPIC
# print('Data transposed')
# data.to_csv('D:/Aaron/Bio/turin/test/gse40279_transposed.csv', header=True)
# data.to_csv('D:/Aaron/Bio/turin/test/gse40279_transposed.txt', sep="\t", header=True)
# data = None

# with open('D:/Aaron/Bio/turin/test/gse40279_transposed.csv', 'a') as f:
#     f.write('\n')
