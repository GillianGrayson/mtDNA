import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, accuracy_score
from sklearn.model_selection import train_test_split

data_path = '../Data/'
data_file_name = 'data_snp_mt.txt'

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

target_pops = ['CHB', 'JPT']    # CHB - Han Chinese in Beijing, China; YRI - Yoruba in Ibadan, Nigeria

reference_pop = 'CHB'   # CHB - Han Chinese in Beijing, China
reference_size = 0.75
reference_list = pop_dict[reference_pop][:int(len(pop_dict[reference_pop])*reference_size)]
reference_frequencies = [0, 0]

# Reference group

f = open(data_path + data_file_name)
header = f.readline().replace('\n', '')
header = header.split(' ')
samples_names = header[15:]

snp_samples = {}
target_samples = []

for sample_name in samples_names:
    if sample_name in reference_list:
        target_samples.append(samples_names.index(sample_name))
        snp_samples[sample_name] = []

line_count = 0
for line in f:
    if line_count % 100 == 0:
        print('Frequencies calculated: ', line_count)

    line = line.replace('\n', '')
    curr_snp_data = line.split(' ')
    snp_chr = curr_snp_data[0]
    snp_gene = curr_snp_data[1]
    snp_pos = curr_snp_data[2]
    snp_name = curr_snp_data[3]
    snp_data = curr_snp_data[15:]

    snp_data = list(snp_data[i] for i in target_samples)

    sample_id = 0
    for sample in snp_data:
        if snp_data[sample_id] == '0':
            reference_frequencies[0] += 1
        elif snp_data[sample_id] == '1':
            reference_frequencies[1] += 1
    sample_id += 1
    line_count += 1

f.close()

reference_frequencies = [freq / sum(reference_frequencies) for freq in reference_frequencies]

# Remaining group

f = open(data_path + data_file_name)
header = f.readline().replace('\n', '')
header = header.split(' ')
samples_names = header[15:]

snp_samples = {}
target_samples = []
genes = []

for sample_name in samples_names:
    for target_pop in target_pops:
        if sample_name in pop_dict[target_pop]:
            target_samples.append(samples_names.index(sample_name))
            snp_samples[sample_name] = []

line_count = 0
for line in f:
    if line_count % 100 == 0:
        print('Lines in mtDNA file: ', line_count)

    line = line.replace('\n', '')
    curr_snp_data = line.split(' ')
    snp_chr = curr_snp_data[0]
    snp_gene = curr_snp_data[1]
    snp_pos = curr_snp_data[2]
    snp_name = curr_snp_data[3]
    snp_data = curr_snp_data[15:]

    snp_data = list(snp_data[i] for i in target_samples)

    sample_id = 0
    genes.append(snp_gene + '_' + snp_pos)
    for sample in snp_samples.keys():
        if snp_data[sample_id] == '0':
            snp_samples[sample].append(1 - reference_frequencies[0])
        elif snp_data[sample_id] == '1':
            snp_samples[sample].append(1 - reference_frequencies[1])
        sample_id += 1
    line_count += 1
f.close()

data_samples = [[] for key in snp_samples.keys()]
data_classes = []
column_names = []
id = 0
for key, value in snp_samples.items():

    if id % 100 == 0:
        print('Samples passed:', id)
    data_samples[id].extend(snp_samples[key])
    if id == 0:
        column_names.extend(genes)

    for target_pop in target_pops:
        if key in pop_dict[target_pop]:
            data_classes.append(target_pop)
    id += 1

snp_samples = {}
genes = []
df = pd.DataFrame(data_samples, columns=column_names)
data_samples = []
df['species'] = data_classes
train, test = train_test_split(df, test_size=0.25)
print('Number of observations in the training data:', len(train))
print('Number of observations in the test data:',len(test))
features = df.columns[:-1]
factor = pd.factorize(train['species'])
y = factor[0]
clf = RandomForestClassifier()
clf.fit(train[features], y)
preds = list(clf.predict(test[features]))
preds = [factor[1][i] for i in preds]

feature_importances = pd.DataFrame(clf.feature_importances_,
                                   index = train[features].columns,
                                   columns=['importance']).sort_values('importance',ascending=False)

print(feature_importances[:10])
print(classification_report(test['species'], preds))
print(accuracy_score(test['species'], preds, normalize=True))
