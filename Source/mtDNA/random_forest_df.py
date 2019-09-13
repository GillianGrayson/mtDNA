import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, accuracy_score
from sklearn.model_selection import train_test_split

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

target_pops = ['CHB', 'YRI']  # CHB - Han Chinese in Beijing, China; YRI - Yoruba in Ibadan, Nigeria

reference_pop = 'CHB'  # CHB - Han Chinese in Beijing, China
reference_size = 0.75
reference_list = pop_dict[reference_pop][:int(len(pop_dict[reference_pop]) * reference_size)]
reference_frequencies = [0, 0, 0, 0, 0, 0]

# Reference group

f = open(data_path + data_file_name_mt)
header = f.readline().replace('\n', '')
header = header.split(' ')
samples_names = header[15:]

target_samples_ids_mt = []
target_samples_names = []
snp_samples_mt = {}

for sample_name in samples_names:
    if sample_name in reference_list:
        target_samples_ids_mt.append(samples_names.index(sample_name))
        target_samples_names.append(sample_name)
        snp_samples_mt[sample_name] = []

for line in f:
    line = line.replace('\n', '')
    curr_snp_data = line.split(' ')
    snp_chr = curr_snp_data[0]
    snp_gene = curr_snp_data[1]
    snp_pos = curr_snp_data[2]
    snp_name = curr_snp_data[3]
    snp_data = curr_snp_data[15:]

    snp_data = list(snp_data[i] for i in target_samples_ids_mt)

    if snp_chr == 'chrMT':
        for id in range(0, len(snp_data)):
            sample_name = samples_names[target_samples_ids_mt[id]]
            snp_samples_mt[sample_name].append(snp_data[id])
f.close()

f = open(data_path + data_file_name)
header = f.readline().replace('\n', '')
header = header.split(' ')
samples_names = header[15:]
target_samples = []
for sample_name in samples_names:
    if sample_name in reference_list and sample_name in snp_samples_mt:
        target_samples.append(samples_names.index(sample_name))

line_count_mt = 0
for line in f:
    if line_count_mt % 100 == 0:
        print('Frequencies calculated: ', line_count_mt)

    line = line.replace('\n', '')
    curr_snp_data = line.split(' ')
    snp_chr = curr_snp_data[0]
    snp_gene = curr_snp_data[1]
    snp_pos = curr_snp_data[2]
    snp_name = curr_snp_data[3]
    snp_data = curr_snp_data[15:]

    snp_data = list(snp_data[i] for i in target_samples)

    sample_id = 0
    for sample in target_samples_names:
        for id in range(0, len(snp_samples_mt[sample])):
            if snp_samples_mt[sample][id] == '0':
                if snp_data[sample_id] == '0|0':
                    reference_frequencies[0] += 1
                elif snp_data[sample_id] == '0|1' or snp_data[sample_id] == '1|0':
                    reference_frequencies[1] += 1
                elif snp_data[sample_id] == '1|1':
                    reference_frequencies[2] += 1
            elif snp_samples_mt[sample][id] == '1':
                if snp_data[sample_id] == '0|0':
                    reference_frequencies[3] += 1
                elif snp_data[sample_id] == '0|1' or snp_data[sample_id] == '1|0':
                    reference_frequencies[4] += 1
                elif snp_data[sample_id] == '1|1':
                    reference_frequencies[5] += 1
        sample_id += 1
    line_count_mt += 1

f.close()

reference_frequencies = [freq / sum(reference_frequencies) for freq in reference_frequencies]
snp_samples_mt = None

# Remaining group

f_mt = open(data_path + data_file_name_mt)
header_mt = f_mt.readline().replace('\n', '')
header_mt = header_mt.split(' ')
samples_names_mt = header_mt[15:]

target_samples_ids_mt = []
target_samples_names_mt = []

for sample_name in samples_names_mt:
    for target_pop in target_pops:
        if sample_name in pop_dict[target_pop]:
            target_samples_names_mt.append(sample_name)
            target_samples_ids_mt.append(samples_names_mt.index(sample_name))

df_ref_mt = []
names_mt = []

line_count_mt = 0
for line in f_mt:
    if line_count_mt % 1000 == 0:
        print('Lines in mtDNA file: ', line_count_mt)

    line = line.replace('\n', '')
    curr_snp_data_mt = line.split(' ')
    snp_chr_mt = curr_snp_data_mt[0]
    snp_gene_mt = curr_snp_data_mt[1]
    snp_pos_mt = curr_snp_data_mt[2]
    snp_name_mt = curr_snp_data_mt[3]
    snp_data_mt = curr_snp_data_mt[15:]

    snp_data_mt = list(snp_data_mt[i] for i in target_samples_ids_mt)

    if snp_chr_mt == 'chrMT':

        name_mt = snp_gene_mt + '_' + snp_pos_mt
        if name_mt not in names_mt:
            names_mt.append(name_mt)

        df_ref_mt.append(snp_data_mt)

    line_count_mt += 1
f_mt.close()


df_ref = pd.DataFrame()

f_nuc = open(data_path + data_file_name)
header_nuc = f_nuc.readline().replace('\n', '')
header_nuc = header_nuc.split(' ')
samples_names_nuc = header_nuc[15:]
target_samples_ids_nuc = []
target_samples_names_nuc = []
for sample_name in samples_names_nuc:
    for target_pop in target_pops:
        if sample_name in pop_dict[target_pop] and sample_name in target_samples_names_mt:
            target_samples_names_nuc.append(sample_name)
            target_samples_ids_nuc.append(samples_names_nuc.index(sample_name))

line_count_nuc = 0
for line in f_nuc:

    if line_count_nuc % 10 == 0:
        print('Lines in nDNA file: ', line_count_nuc)

    line = line.replace('\n', '')
    curr_snp_data_nuc = line.split(' ')
    snp_chr_nuc = curr_snp_data_nuc[0]
    snp_gene_nuc = curr_snp_data_nuc[1]
    snp_pos_nuc = curr_snp_data_nuc[2]
    snp_name_nuc = curr_snp_data_nuc[3]
    snp_data_nuc = curr_snp_data_nuc[15:]

    snp_data_nuc = list(snp_data_nuc[i] for i in target_samples_ids_nuc)

    for mt_id in range(0, len(names_mt)):

        name_mt = names_mt[mt_id]

        combination_name = name_mt + '_' + snp_gene_nuc + '_' + snp_name_nuc
        combination_data = []

        for id in range(0, len(snp_data_nuc)):
            if df_ref_mt[mt_id][id] == '0':
                if snp_data_nuc[id] == '0|0':
                    combination_data.append(1 - reference_frequencies[0])
                elif snp_data_nuc[id] == '0|1' or snp_data_nuc[id] == '1|0':
                    combination_data.append(1 - reference_frequencies[1])
                elif snp_data_nuc[id] == '1|1':
                    combination_data.append(1 - reference_frequencies[2])
            elif df_ref_mt[mt_id][id] == '1':
                if snp_data_nuc[id] == '0|0':
                    combination_data.append(1 - reference_frequencies[3])
                elif snp_data_nuc[id] == '0|1' or snp_data_nuc[id] == '1|0':
                    combination_data.append(1 - reference_frequencies[4])
                elif snp_data_nuc[id] == '1|1':
                    combination_data.append(1 - reference_frequencies[5])

        df_ref[combination_name] = combination_data

    line_count_nuc += 1

f_nuc.close()

data_classes = []
for item in target_samples_names_mt:
    for target_pop in target_pops:
        if item in pop_dict[target_pop]:
            data_classes.append(target_pop)

df_ref['species'] = data_classes
train, test = train_test_split(df_ref, test_size=0.25)
print('Number of observations in the training data:', len(train))
print('Number of observations in the test data:', len(test))
features = df_ref.columns[:-1]
factor = pd.factorize(train['species'])
y = factor[0]
clf = RandomForestClassifier()
clf.fit(train[features], y)
preds = list(clf.predict(test[features]))
preds = [factor[1][i] for i in preds]

feature_importances = pd.DataFrame(clf.feature_importances_,
                                   index=train[features].columns,
                                   columns=['importance']).sort_values('importance', ascending=False)

print(feature_importances[:10])
print(classification_report(test['species'], preds))
print(accuracy_score(test['species'], preds, normalize=True))
