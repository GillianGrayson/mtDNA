import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report

data_path = '../Data/'
data_file_name = 'data_snp_test_genes_short.txt'
data_file_name_mt = 'data_snp_mt.txt'

pop_file_name = 's_pop.txt'
pop_dict = {}
f = open(data_path + pop_file_name)
f.readline()
for line in f:
    line = line.replace('\n', '')
    curr_pop_data = line.split('\t')
    if curr_pop_data[1] in pop_dict:
        pop_dict[curr_pop_data[1]].append(curr_pop_data[0])
    else:
        pop_dict[curr_pop_data[1]] = [curr_pop_data[0]]
f.close()

target_pops = ['CHB', 'YRI'] # CHB - Han Chinese in Beijing, China; YRI - Yoruba in Ibadan, Nigeria

f = open(data_path + data_file_name_mt)
header = f.readline().replace('\n', '')
header = header.split(' ')
samples_names = header[15:]

snp_samples = {}
target_samples_mt = []
snp_samples_mt = {}
genes_mt = {}
genes_combinations = {}

for sample_name in samples_names:
    for target_pop in target_pops:
        if sample_name in pop_dict[target_pop]:
            target_samples_mt.append(samples_names.index(sample_name))
            snp_samples[sample_name] = []
            snp_samples_mt[sample_name] = None
            genes_mt[sample_name] = []
            genes_combinations[sample_name] = []


for line in f:
    line = line.replace('\n', '')
    curr_snp_data = line.split(' ')
    snp_chr = curr_snp_data[0]
    snp_gene = curr_snp_data[1]
    snp_pos = curr_snp_data[2]
    snp_name = curr_snp_data[3]
    snp_data = curr_snp_data[15:]

    snp_data = list(snp_data[i] for i in target_samples_mt)

    if snp_chr == 'chrMT':
        for id in range(0, len(snp_data)):
            sample_name = samples_names[target_samples_mt[id]]
            snp_samples_mt[sample_name] = snp_data[id]
            genes_mt[sample_name].append(snp_gene + '_' + snp_name)
f.close()

f = open(data_path + data_file_name)
header = f.readline().replace('\n', '')
header = header.split(' ')
samples_names = header[15:]
target_samples = []
for sample_name in samples_names:
    for target_pop in target_pops:
        if sample_name in pop_dict[target_pop] and sample_name in snp_samples_mt:
            target_samples.append(samples_names.index(sample_name))
for line in f:
    line = line.replace('\n', '')
    curr_snp_data = line.split(' ')
    snp_chr = curr_snp_data[0]
    snp_gene = curr_snp_data[1]
    snp_pos = curr_snp_data[2]
    snp_name = curr_snp_data[3]
    snp_data = curr_snp_data[15:]

    snp_data = list(snp_data[i] for i in target_samples)

    for sample in snp_samples.keys():
        for id in range(0, len(snp_samples_mt)):
            sample_name = samples_names[target_samples[id]]
            genes_combinations[sample_name].append(
                [genes_mt[sample_name][id] + '_' + snp_gene + '_' + snp_name + '_0_0|0',
                 genes_mt[sample_name][id] + '_' + snp_gene + '_' + snp_name + '_0_1|0',
                 genes_mt[sample_name][id] + '_' + snp_gene + '_' + snp_name + '_0_1|1',
                 genes_mt[sample_name][id] + '_' + snp_gene + '_' + snp_name + '_1_0|0',
                 genes_mt[sample_name][id] + '_' + snp_gene + '_' + snp_name + '_1_1|0',
                 genes_mt[sample_name][id] + '_' + snp_gene + '_' + snp_name + '_1_1|1'])
            if snp_samples_mt[sample_name] == '0':
                if snp_data[id] == '0|0':
                    snp_samples[sample_name].append([1, 0, 0, 0, 0, 0])  # 0 and 0|0 is [1, 0, 0, 0, 0, 0]
                elif snp_data[id] == '0|1' or snp_data[id] == '1|0':
                    snp_samples[sample_name].append([0, 1, 0, 0, 0, 0])  # 0 and 0|1 or 1|0 is [0, 1, 0, 0, 0, 0]
                elif snp_data[id] == '1|1':
                    snp_samples[sample_name].append([0, 0, 1, 0, 0, 0])  # 0 and 1|1 is [0, 0, 1, 0, 0, 0]
            elif snp_samples_mt[sample_name] == '1':
                if snp_data[id] == '0|0':
                    snp_samples[sample_name].append([0, 0, 0, 1, 0, 0])  # 1 and 0|0 is [0, 0, 0, 1, 0, 0]
                elif snp_data[id] == '0|1' or snp_data[id] == '1|0':
                    snp_samples[sample_name].append([0, 0, 0, 0, 1, 0])  # 1 and 0|1 or 1|0 is [0, 0, 0, 0, 1, 0]
                elif snp_data[id] == '1|1':
                    snp_samples[sample_name].append([0, 0, 0, 0, 0, 1])  # 1 and 1|1 is [0, 0, 0, 0, 0, 1]
f.close()

data_samples = []
data_classes = []
column_names = []
id = 0
for key, value in snp_samples.items():
    data_samples.append(snp_samples[key][0])

    if id == 0:
        column_names.extend(genes_combinations[key][0])

    for target_pop in target_pops:
        if key in pop_dict[target_pop]:
            data_classes.append(target_pop)

    for i in range(1, len(value)):
        data_samples[id].extend(value[i])
        if id == 0:
            column_names.extend(genes_combinations[key][i])
    id += 1

snp_samples = {}
genes_combinations = {}
df = pd.DataFrame(data_samples, columns=column_names)
data_samples = []
df['species'] = data_classes
df['is_train'] = np.random.uniform(0, 1, len(df)) <= .75
train, test = df[df['is_train']==True], df[df['is_train']==False]
print('Number of observations in the training data:', len(train))
print('Number of observations in the test data:',len(test))
features = df.columns[:-2]
y = pd.factorize(train['species'])[0]
clf = RandomForestClassifier(n_jobs=2, random_state=0)
clf.fit(train[features], y)
preds = [target_pops[i] for i in list(clf.predict(test[features]))]
print(classification_report(test['species'], preds))

ololo = 1