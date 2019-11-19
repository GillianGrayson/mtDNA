import os
import itertools
import numpy as np
import pandas as pd
from tqdm import tqdm
from collections import Counter
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_validate

data_path = 'D:/Aaron/Bio/tibet/Data/'
result_path = 'D:/Aaron/Bio/tibet/Result/'

raw_data = []
data_classes = []
for filename in os.listdir(data_path):
    f = open(data_path + filename, 'r')
    raw_data.append([line.rstrip() for line in f][1::2])
    data_classes.append(filename[:-6])
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
classes = []
for i in range(0, len(raw_data)):
    classes += [data_classes[i], ] * len(raw_data[i])

factor = pd.factorize(classes)
y = factor[0]
clf = RandomForestClassifier(n_estimators=500)
output = cross_validate(clf, df_main, y, cv=10, scoring='accuracy', return_estimator=True)
accuracy = np.mean(output['test_score'])
print('8-class Classification Accuracy: ' + str(accuracy))

result_file_name = 'classification.txt'
if not os.path.exists(result_path):
    os.makedirs(result_path)
f = open(result_path + result_file_name, 'w')
f.write('8-class Classification Accuracy: ' + str(accuracy) + '\n')

for subset in itertools.combinations(raw_data, 2):
    test_data = subset
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
    classes = [data_classes[raw_data.index(test_data[0])], ] * len(test_data[0]) + \
              [data_classes[raw_data.index(test_data[1])], ] * len(test_data[1])

    factor = pd.factorize(classes)
    y = factor[0]
    clf = RandomForestClassifier(n_estimators=500)
    output = cross_validate(clf, df_main, y, cv=10, scoring='accuracy', return_estimator=True)
    accuracy = np.mean(output['test_score'])
    print(data_classes[raw_data.index(test_data[0])] + ' vs ' +  data_classes[raw_data.index(test_data[1])] +
          ' Binary Classification Accuracy: ' + str(accuracy))
    f.write(data_classes[raw_data.index(test_data[0])] + ' vs ' +  data_classes[raw_data.index(test_data[1])] +
          ' Binary Classification Accuracy: ' + str(accuracy) + '\n')
f.close()
