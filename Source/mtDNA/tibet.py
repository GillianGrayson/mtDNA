import os
import numpy as np
from tqdm import tqdm
from collections import Counter

data_path = 'D:/Aaron/Bio/tibet/Data/'

raw_data = []
for filename in os.listdir(data_path):
    f = open(data_path + filename,'r')
    raw_data.append([line.rstrip() for line in f][1::2])
    f.close()

num_nuc = len(raw_data[0][0])
num_persons = int(np.sum([len(raw_data[i]) for i in range(0, len(raw_data))]))
df_main = np.empty(shape=(num_persons, num_nuc), dtype=np.int)

for nuc_id in tqdm(range(0, num_nuc)):
    curr_nuc = [raw_data[i][j][nuc_id] for i in range(0, len(raw_data)) for j in range(0, len(raw_data[i]))]
    count_dict = Counter(curr_nuc)
    count_dict.most_common()
    if len(count_dict) == 1:
        for person_id in range(0, num_persons):
            df_main[person_id, nuc_id] = 0
    else:
        for person_id in range(0, num_persons):
            df_main[person_id, nuc_id] = list(count_dict.keys()).index(curr_nuc[person_id])

df_main = df_main[:, ~np.all(df_main[1:] == df_main[:-1], axis=0)]
olo =0
