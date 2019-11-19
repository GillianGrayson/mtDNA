import os
import numpy as np
import pandas as pd
from tqdm import tqdm
from collections import Counter
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_validate

data_path = 'D:/Aaron/Bio/tibet/Data/'

raw_data = []
for filename in os.listdir(data_path):
    f = open(data_path + filename, 'r')
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
classes = ['0-500', ] * len(raw_data[0]) + ['501-1000', ] * len(raw_data[1]) + \
          ['1001-1500', ] * len(raw_data[2]) + ['1501-2000', ] * len(raw_data[3]) + \
          ['2001-2500', ] * len(raw_data[4]) + ['2501-3000', ] * len(raw_data[5]) + \
          ['3001-4000', ] * len(raw_data[6]) + ['4001', ] * len(raw_data[7])

factor = pd.factorize(classes)
y = factor[0]
clf = RandomForestClassifier(n_estimators=500)
output = cross_validate(clf, df_main, y, cv=10, scoring='accuracy', return_estimator=True)
accuracy = np.mean(output['test_score'])
print('Classification Accuracy: ' + str(accuracy))

test_data = [raw_data[0], raw_data[-1]]
num_persons = int(np.sum([len(test_data[i]) for i in range(0, len(test_data))]))
df_main = np.empty(shape=(num_persons, num_nuc), dtype=np.int)
for nuc_id in tqdm(range(0, num_nuc)):
    curr_nuc = [test_data[i][j][nuc_id] for i in range(0, len(test_data)) for j in range(0, len(test_data[i]))]
    count_dict = Counter(curr_nuc)
    count_dict.most_common()
    if len(count_dict) == 1:
        for person_id in range(0, num_persons):
            df_main[person_id, nuc_id] = 0
    else:
        for person_id in range(0, num_persons):
            df_main[person_id, nuc_id] = list(count_dict.keys()).index(curr_nuc[person_id])

df_main = df_main[:, ~np.all(df_main[1:] == df_main[:-1], axis=0)]
classes = ['0-500', ] * len(test_data[0]) + ['4001', ] * len(test_data[1])

factor = pd.factorize(classes)
y = factor[0]
clf = RandomForestClassifier(n_estimators=500)
output = cross_validate(clf, df_main, y, cv=10, scoring='accuracy', return_estimator=True)
accuracy = np.mean(output['test_score'])
print('Binary Classification Accuracy: ' + str(accuracy))
