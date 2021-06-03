from Source.mtDNA.tibet.functions.file_system import get_path
import pandas as pd
import numpy as np
from tqdm import tqdm
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_validate
from collections import Counter


path = get_path()
alignment_data_path = path + '/Data/alignment/low_data/'
info_data_path = path + '/Data/alignment/info/'
result_data_path = path + '/Result/align/low_data/'

# Read data file
data_dict = {}
f = open(alignment_data_path + 'all_low.txt', 'r')
for line in f:
    if line.startswith('>'):
        curr_key = line.rstrip().split(' ')[0][1:]
        data_dict[curr_key] = ''
    elif prev_line.startswith('>'):
        data_dict[curr_key] = line.rstrip()
    else:
        data_dict[curr_key] += line.rstrip()
    prev_line = line
f.close()

rcrs = data_dict['rcrs']
data_dict.pop('rcrs', None)

# Read subjects groups file
subjects_data = pd.read_excel(info_data_path + 'subjects_all.xlsx').to_dict('list')
high_classes = ['Andes', 'Tibetan', 'Ethiopia', '4001']
subject_info = {high_class: [] for high_class in high_classes}
for subject in data_dict:
    if subject.split('_')[0] in high_classes:
        subject_info[subject.split('_')[0]].append(subject)
        subject_info['4001'].append(subject)
    else:
        subject_id = subjects_data['subject'].index(subject)
        subject_group = subjects_data['height'][subject_id]
        if subject_group not in subject_info:
            subject_info[subject_group] = [subject]
        else:
            subject_info[subject_group].append(subject)

# Align data to rCRS
num_nucleotides = len(rcrs.replace('~', ''))
num_subjects = len(data_dict)
table = np.empty(shape=(num_subjects, num_nucleotides), dtype=object)
num_gaps = 0
for nucleotide_id in tqdm(range(0, len(rcrs))):
    curr_nucleotide = rcrs[nucleotide_id]
    if curr_nucleotide == '~':
        num_gaps += 1
        for subject_id in range(0, num_subjects):
            curr_subject = list(data_dict.keys())[subject_id]
            if data_dict[curr_subject][nucleotide_id] != '~':
                table[subject_id, nucleotide_id - num_gaps] = table[subject_id, nucleotide_id - num_gaps] + \
                                                              data_dict[curr_subject][nucleotide_id]
    else:
        for subject_id in range(0, num_subjects):
            curr_subject = list(data_dict.keys())[subject_id]
            table[subject_id, nucleotide_id - num_gaps] = data_dict[curr_subject][nucleotide_id]
positions = list(range(0, num_nucleotides))
rcrs = rcrs.replace('~', '')

# Read phylotrees data
phylotrees = pd.read_excel(info_data_path + 'phylotrees.xlsx').to_dict('list')
phylotrees_positions = []
for position in phylotrees['position']:
    if isinstance(position, int):
        position = int(position) - 1
        if position not in phylotrees_positions:
            phylotrees_positions.append(position)
    else:
        start = int(position.split('-')[0])
        finish = int(position.split('-')[1])
        positions = range(start, finish + 1)
        for curr_position in positions:
            curr_position -= 1
            if curr_position not in phylotrees_positions:
                phylotrees_positions.append(curr_position)
phylotrees_positions.sort()

# Aligned data with phylotrees
table_phylo = np.delete(table, phylotrees_positions, 1)
index_phylotrees = [i for j, i in enumerate(list(range(0, num_nucleotides))) if j not in phylotrees_positions]

# Encode data as 0/1
table_code = np.empty(shape=(num_subjects, len(index_phylotrees)), dtype=int)
for nucleotide_id in tqdm(range(0, len(index_phylotrees))):
    rcrs_id = index_phylotrees[nucleotide_id]
    curr_rcrs_nucleotide = rcrs[rcrs_id]
    for subject_id in range(0, num_subjects):
        if table_phylo[subject_id, nucleotide_id] == curr_rcrs_nucleotide:
            table_code[subject_id, nucleotide_id] = 0
        else:
            table_code[subject_id, nucleotide_id] = 1
"""
# Encode data using frequency
table_code = np.empty(shape=(num_subjects, len(index_phylotrees)), dtype=int)
for nucleotide_id in tqdm(range(0, len(index_phylotrees))):
    rcrs_id = index_phylotrees[nucleotide_id]
    curr_rcrs_nucleotide = rcrs[rcrs_id]
    curr_nucleotide = table_phylo[:, nucleotide_id]
    count_dict = Counter(curr_nucleotide).most_common()
    count_dict = dict(count_dict)
    variants_list = list(count_dict.keys())
    for subject_id in range(0, num_subjects):
        curr_subject_nucleotide = table_phylo[subject_id, nucleotide_id]
        table_code[subject_id, nucleotide_id] = variants_list.index(curr_subject_nucleotide)
"""
all_valuable_features = []
# Low-High classification
curr_exp_classes = ['0-500', '4001']
curr_exp_indexes = []
classes = []
for low_high_class in curr_exp_classes:
    for subject in subject_info[low_high_class]:
        curr_exp_indexes.append(list(data_dict.keys()).index(subject))
        classes.append(low_high_class)
factor = pd.factorize(classes)
y = factor[0]
curr_table = table_code[curr_exp_indexes, :]
rf_classifier = RandomForestClassifier(n_estimators=500)
output = cross_validate(rf_classifier, curr_table, y, cv=5, scoring='accuracy', return_estimator=True)
mean_accuracy = np.mean(output['test_score'])
features_dict = dict((key, []) for key in index_phylotrees)
for idx, estimator in enumerate(output['estimator']):
    feature_importances = pd.DataFrame(estimator.feature_importances_,
                                       index=index_phylotrees,
                                       columns=['importance']).sort_values('importance', ascending=False)
    features_names = list(feature_importances.index.values)
    features_values = list(feature_importances.values)
    for i in range(0, len(features_names)):
        features_dict[features_names[i]].append(features_values[i][0])
for key in features_dict.keys():
    features_dict[key] = np.mean(features_dict[key])
features_dict = {k: v for k, v in sorted(features_dict.items(), reverse=True, key=lambda x: x[1])}
features_rating = []
for feature in features_dict:
    if features_dict[feature] > 0:
        features_rating.append(feature)
all_valuable_features.extend(features_rating)

# Save Low-High classification results
f = open(result_data_path + '0-500_4001.txt', 'w')
f.write(str(mean_accuracy) + '\n')
for feature in features_rating:
    f.write(str(feature + 1) + '\n')
f.close()

# Low-Andes classification
curr_exp_classes = ['0-500', 'Andes']
curr_exp_indexes = []
classes = []
for low_high_class in curr_exp_classes:
    for subject in subject_info[low_high_class]:
        curr_exp_indexes.append(list(data_dict.keys()).index(subject))
        classes.append(low_high_class)
factor = pd.factorize(classes)
y = factor[0]
curr_table = table_code[curr_exp_indexes, :]
rf_classifier = RandomForestClassifier(n_estimators=500)
output = cross_validate(rf_classifier, curr_table, y, cv=5, scoring='accuracy', return_estimator=True)
mean_accuracy = np.mean(output['test_score'])
features_dict = dict((key, []) for key in index_phylotrees)
for idx, estimator in enumerate(output['estimator']):
    feature_importances = pd.DataFrame(estimator.feature_importances_,
                                       index=index_phylotrees,
                                       columns=['importance']).sort_values('importance', ascending=False)
    features_names = list(feature_importances.index.values)
    features_values = list(feature_importances.values)
    for i in range(0, len(features_names)):
        features_dict[features_names[i]].append(features_values[i][0])
for key in features_dict.keys():
    features_dict[key] = np.mean(features_dict[key])
features_dict = {k: v for k, v in sorted(features_dict.items(), reverse=True, key=lambda x: x[1])}
features_rating = []
for feature in features_dict:
    if features_dict[feature] > 0:
        features_rating.append(feature)

# Save Low-Andes classification results
f = open(result_data_path + '0-500_Andes.txt', 'w')
f.write(str(mean_accuracy) + '\n')
for feature in features_rating:
    f.write(str(feature + 1) + '\n')
f.close()
all_valuable_features.extend(features_rating)

# Low-Tibetan classification
curr_exp_classes = ['0-500', 'Tibetan']
curr_exp_indexes = []
classes = []
for low_high_class in curr_exp_classes:
    for subject in subject_info[low_high_class]:
        curr_exp_indexes.append(list(data_dict.keys()).index(subject))
        classes.append(low_high_class)
factor = pd.factorize(classes)
y = factor[0]
curr_table = table_code[curr_exp_indexes, :]
rf_classifier = RandomForestClassifier(n_estimators=500)
output = cross_validate(rf_classifier, curr_table, y, cv=5, scoring='accuracy', return_estimator=True)
mean_accuracy = np.mean(output['test_score'])
features_dict = dict((key, []) for key in index_phylotrees)
for idx, estimator in enumerate(output['estimator']):
    feature_importances = pd.DataFrame(estimator.feature_importances_,
                                       index=index_phylotrees,
                                       columns=['importance']).sort_values('importance', ascending=False)
    features_names = list(feature_importances.index.values)
    features_values = list(feature_importances.values)
    for i in range(0, len(features_names)):
        features_dict[features_names[i]].append(features_values[i][0])
for key in features_dict.keys():
    features_dict[key] = np.mean(features_dict[key])
features_dict = {k: v for k, v in sorted(features_dict.items(), reverse=True, key=lambda x: x[1])}
features_rating = []
for feature in features_dict:
    if features_dict[feature] > 0:
        features_rating.append(feature)

# Save Low-Tibetan classification results
f = open(result_data_path + '0-500_Tibetan.txt', 'w')
f.write(str(mean_accuracy) + '\n')
for feature in features_rating:
    f.write(str(feature + 1) + '\n')
f.close()
all_valuable_features.extend(features_rating)

# Low-Ethiopia classification
curr_exp_classes = ['0-500', 'Ethiopia']
curr_exp_indexes = []
classes = []
for low_high_class in curr_exp_classes:
    for subject in subject_info[low_high_class]:
        curr_exp_indexes.append(list(data_dict.keys()).index(subject))
        classes.append(low_high_class)
factor = pd.factorize(classes)
y = factor[0]
curr_table = table_code[curr_exp_indexes, :]
rf_classifier = RandomForestClassifier(n_estimators=500)
output = cross_validate(rf_classifier, curr_table, y, cv=5, scoring='accuracy', return_estimator=True)
mean_accuracy = np.mean(output['test_score'])
features_dict = dict((key, []) for key in index_phylotrees)
for idx, estimator in enumerate(output['estimator']):
    feature_importances = pd.DataFrame(estimator.feature_importances_,
                                       index=index_phylotrees,
                                       columns=['importance']).sort_values('importance', ascending=False)
    features_names = list(feature_importances.index.values)
    features_values = list(feature_importances.values)
    for i in range(0, len(features_names)):
        features_dict[features_names[i]].append(features_values[i][0])
for key in features_dict.keys():
    features_dict[key] = np.mean(features_dict[key])
features_dict = {k: v for k, v in sorted(features_dict.items(), reverse=True, key=lambda x: x[1])}
features_rating = []
for feature in features_dict:
    if features_dict[feature] > 0:
        features_rating.append(feature)

# Save Low-Ethiopia classification results
f = open(result_data_path + '0-500_Ethiopia.txt', 'w')
f.write(str(mean_accuracy) + '\n')
for feature in features_rating:
    f.write(str(feature + 1) + '\n')
f.close()
all_valuable_features.extend(features_rating)

# Tibetan-Andes classification
curr_exp_classes = ['Tibetan', 'Andes']
curr_exp_indexes = []
classes = []
for low_high_class in curr_exp_classes:
    for subject in subject_info[low_high_class]:
        curr_exp_indexes.append(list(data_dict.keys()).index(subject))
        classes.append(low_high_class)
factor = pd.factorize(classes)
y = factor[0]
curr_table = table_code[curr_exp_indexes, :]
rf_classifier = RandomForestClassifier(n_estimators=500)
output = cross_validate(rf_classifier, curr_table, y, cv=5, scoring='accuracy', return_estimator=True)
mean_accuracy = np.mean(output['test_score'])
features_dict = dict((key, []) for key in index_phylotrees)
for idx, estimator in enumerate(output['estimator']):
    feature_importances = pd.DataFrame(estimator.feature_importances_,
                                       index=index_phylotrees,
                                       columns=['importance']).sort_values('importance', ascending=False)
    features_names = list(feature_importances.index.values)
    features_values = list(feature_importances.values)
    for i in range(0, len(features_names)):
        features_dict[features_names[i]].append(features_values[i][0])
for key in features_dict.keys():
    features_dict[key] = np.mean(features_dict[key])
features_dict = {k: v for k, v in sorted(features_dict.items(), reverse=True, key=lambda x: x[1])}
features_rating = []
for feature in features_dict:
    if features_dict[feature] > 0:
        features_rating.append(feature)

# Save Tibetan-Andes classification results
f = open(result_data_path + 'Tibetan_Andes.txt', 'w')
f.write(str(mean_accuracy) + '\n')
for feature in features_rating:
    f.write(str(feature + 1) + '\n')
f.close()
all_valuable_features.extend(features_rating)

# Tibetan-Ethiopia classification
curr_exp_classes = ['Tibetan', 'Ethiopia']
curr_exp_indexes = []
classes = []
for low_high_class in curr_exp_classes:
    for subject in subject_info[low_high_class]:
        curr_exp_indexes.append(list(data_dict.keys()).index(subject))
        classes.append(low_high_class)
factor = pd.factorize(classes)
y = factor[0]
curr_table = table_code[curr_exp_indexes, :]
rf_classifier = RandomForestClassifier(n_estimators=500)
output = cross_validate(rf_classifier, curr_table, y, cv=5, scoring='accuracy', return_estimator=True)
mean_accuracy = np.mean(output['test_score'])
features_dict = dict((key, []) for key in index_phylotrees)
for idx, estimator in enumerate(output['estimator']):
    feature_importances = pd.DataFrame(estimator.feature_importances_,
                                       index=index_phylotrees,
                                       columns=['importance']).sort_values('importance', ascending=False)
    features_names = list(feature_importances.index.values)
    features_values = list(feature_importances.values)
    for i in range(0, len(features_names)):
        features_dict[features_names[i]].append(features_values[i][0])
for key in features_dict.keys():
    features_dict[key] = np.mean(features_dict[key])
features_dict = {k: v for k, v in sorted(features_dict.items(), reverse=True, key=lambda x: x[1])}
features_rating = []
for feature in features_dict:
    if features_dict[feature] > 0:
        features_rating.append(feature)

# Save Tibetan-Ethiopia classification results
f = open(result_data_path + 'Tibetan_Ethiopia.txt', 'w')
f.write(str(mean_accuracy) + '\n')
for feature in features_rating:
    f.write(str(feature + 1) + '\n')
f.close()
all_valuable_features.extend(features_rating)

# Andes-Ethiopia classification
curr_exp_classes = ['Andes', 'Ethiopia']
curr_exp_indexes = []
classes = []
for low_high_class in curr_exp_classes:
    for subject in subject_info[low_high_class]:
        curr_exp_indexes.append(list(data_dict.keys()).index(subject))
        classes.append(low_high_class)
factor = pd.factorize(classes)
y = factor[0]
curr_table = table_code[curr_exp_indexes, :]
rf_classifier = RandomForestClassifier(n_estimators=500)
output = cross_validate(rf_classifier, curr_table, y, cv=5, scoring='accuracy', return_estimator=True)
mean_accuracy = np.mean(output['test_score'])
features_dict = dict((key, []) for key in index_phylotrees)
for idx, estimator in enumerate(output['estimator']):
    feature_importances = pd.DataFrame(estimator.feature_importances_,
                                       index=index_phylotrees,
                                       columns=['importance']).sort_values('importance', ascending=False)
    features_names = list(feature_importances.index.values)
    features_values = list(feature_importances.values)
    for i in range(0, len(features_names)):
        features_dict[features_names[i]].append(features_values[i][0])
for key in features_dict.keys():
    features_dict[key] = np.mean(features_dict[key])
features_dict = {k: v for k, v in sorted(features_dict.items(), reverse=True, key=lambda x: x[1])}
features_rating = []
for feature in features_dict:
    if features_dict[feature] > 0:
        features_rating.append(feature)

# Save Andes-Ethiopia classification results
f = open(result_data_path + 'Andes_Ethiopia.txt', 'w')
f.write(str(mean_accuracy) + '\n')
for feature in features_rating:
    f.write(str(feature + 1) + '\n')
f.close()
all_valuable_features.extend(features_rating)

all_valuable_features = list(set(all_valuable_features))
all_valuable_features.sort()

# TEA classification
curr_exp_classes = ['Tibetan', 'Ethiopia',  'Andes']
curr_exp_indexes = []
classes = []
for low_high_class in curr_exp_classes:
    for subject in subject_info[low_high_class]:
        curr_exp_indexes.append(list(data_dict.keys()).index(subject))
        classes.append(low_high_class)
factor = pd.factorize(classes)
y = factor[0]
curr_table = table_code[curr_exp_indexes, :]
rf_classifier = RandomForestClassifier(n_estimators=500)
output = cross_validate(rf_classifier, curr_table, y, cv=5, scoring='accuracy', return_estimator=True)
mean_accuracy = np.mean(output['test_score'])
features_dict = dict((key, []) for key in index_phylotrees)
for idx, estimator in enumerate(output['estimator']):
    feature_importances = pd.DataFrame(estimator.feature_importances_,
                                       index=index_phylotrees,
                                       columns=['importance']).sort_values('importance', ascending=False)
    features_names = list(feature_importances.index.values)
    features_values = list(feature_importances.values)
    for i in range(0, len(features_names)):
        features_dict[features_names[i]].append(features_values[i][0])
for key in features_dict.keys():
    features_dict[key] = np.mean(features_dict[key])
features_dict = {k: v for k, v in sorted(features_dict.items(), reverse=True, key=lambda x: x[1])}
features_rating = []
for feature in features_dict:
    if features_dict[feature] > 0:
        features_rating.append(feature)

# Save TEA classification results
f = open(result_data_path + 'Tibetan_Andes_Ethiopia.txt', 'w')
f.write(str(mean_accuracy) + '\n')
for feature in features_rating:
    f.write(str(feature + 1) + '\n')
f.close()
all_valuable_features.extend(features_rating)

# Calculate mutation statistic
all_classes = ['0-500', 'Tibetan', 'Andes', 'Ethiopia']
all_classes_indexes = {key: [] for key in all_classes}
for curr_class in all_classes:
    for subject in subject_info[curr_class]:
        all_classes_indexes[curr_class].append(list(data_dict.keys()).index(subject))
mutation_stat_dict = {'Position': []}
classes_extended = []
for curr_class in all_classes:
    classes_extended.append(curr_class + ' Main')
    classes_extended.append(curr_class + ' Main Freq')
    classes_extended.append(curr_class + ' Minor')
    classes_extended.append(curr_class + ' Minor Freq')
mutation_stat_dict.update({curr_class: [] for curr_class in classes_extended})
for position in all_valuable_features:
    mutation_stat_dict['Position'].append(position + 1)
    for curr_class in all_classes:
        curr_nucleotide = table[all_classes_indexes[curr_class], position]
        count_dict = Counter(curr_nucleotide).most_common()
        count_dict = dict(count_dict)
        variants_list = list(count_dict.keys())
        mutation_stat_dict[curr_class + ' Main'].append(variants_list[0])
        mutation_stat_dict[curr_class + ' Main Freq'].append(count_dict[variants_list[0]] / len(curr_nucleotide))
        if len(variants_list) == 1:
            mutation_stat_dict[curr_class + ' Minor'].append('')
            mutation_stat_dict[curr_class + ' Minor Freq'].append(0.0)
        if len(variants_list) > 1:
            mutation_stat_dict[curr_class + ' Minor'].append(variants_list[1])
            mutation_stat_dict[curr_class + ' Minor Freq'].append(count_dict[variants_list[1]] / len(curr_nucleotide))

# Save mutation statistic
df = pd.DataFrame(mutation_stat_dict)
df.to_excel(result_data_path + 'stat.xlsx', index=False)
