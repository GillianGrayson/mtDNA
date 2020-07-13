import os
import numpy as np
import pandas as pd
from tqdm import tqdm
from collections import Counter
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_validate


def read_data(data_path):
    raw_data = []
    subjects = []
    data_classes = []
    for filename in os.listdir(data_path):
        if filename.endswith('fasta') or filename.endswith('fa'):
            f = open(data_path + filename, 'r')
            raw_data.append([line.rstrip() for line in f][1::2])
            f = open(data_path + filename, 'r')
            subjects.append([line.rstrip().split(' ')[0][1:] for line in f][0::2])
            if filename.endswith('fasta'):
                data_classes.append(filename[:-6])
            else:
                data_classes.append(filename[:-3])
            f.close()
    return raw_data, subjects, data_classes


def get_region_info(data_path):
    regions = {'region': [], 'start': [], 'finish': []}
    f = open(data_path + 'regions.txt', 'r')
    for line in f:
        line_list = line.rstrip().split('\t')
        regions['region'].append(line_list[0])
        regions['start'].append(int(line_list[1]))
        regions['finish'].append(int(line_list[2]))
    f.close()
    return regions


def subset_subjects(raw_data, initial_classes, subset_classes):
    subset_data = {key: [] for key in subset_classes}
    subject_classes = []
    for key in subset_classes:
        for curr_class in subset_classes[key]:
            class_id = initial_classes.index(curr_class)
            subset_data[key].extend(raw_data[class_id])
        subject_classes.extend([key] * len(subset_data[key]))
    return subset_data, subject_classes


def create_classes_table(raw_data):
    data_classes = list(raw_data.keys())
    num_nucleotides = len(raw_data[data_classes[0]][0])
    num_persons = int(np.sum([len(raw_data[data_class]) for data_class in data_classes]))
    table = np.empty(shape=(num_persons, num_nucleotides), dtype=np.int)
    mut_positions = []

    for nucleotide_id in tqdm(range(0, num_nucleotides)):
        curr_nuc = [raw_data[data_class][person_id][nucleotide_id] for data_class in data_classes for person_id in
                    range(0, len(raw_data[data_class]))]
        count_dict = Counter(curr_nuc).most_common()
        count_dict = dict(count_dict)
        if '-' in count_dict:
            if count_dict['-'] / len(curr_nuc) > 0.9:
                continue
        if len(count_dict) == 1:
            for person_id in range(0, num_persons):
                table[person_id, nucleotide_id] = 0
        else:
            for person_id in range(0, num_persons):
                if list(count_dict.keys()).index(curr_nuc[person_id]) == 0:
                    table[person_id, nucleotide_id] = 0
                else:
                    table[person_id, nucleotide_id] = 1
        if len(np.unique(table[:, nucleotide_id])) > 1:
            mut_positions.append(nucleotide_id)

    table = table[:, ~np.all(table[1:] == table[:-1], axis=0)]
    return table, mut_positions


def run_sequential_random_forest(table, classes, positions, num_runs):
    factor = pd.factorize(classes)
    y = factor[0]
    clf = RandomForestClassifier(n_estimators=500)
    output = cross_validate(clf, table, y, cv=10, scoring='accuracy', return_estimator=True)

    features_dict = dict((key, []) for key in positions)
    for idx, estimator in enumerate(output['estimator']):
        feature_importances = pd.DataFrame(estimator.feature_importances_,
                                           index=positions,
                                           columns=['importance']).sort_values('importance', ascending=False)

        features_names = list(feature_importances.index.values)
        features_values = list(feature_importances.values)
        for i in range(0, len(features_names)):
            features_dict[features_names[i]].append(features_values[i][0])
    for key in features_dict.keys():
        features_dict[key] = np.mean(features_dict[key])
    features_dict = {k: v for k, v in sorted(features_dict.items(), reverse=True, key=lambda x: x[1])}

    features_rating = list(features_dict.keys())
    accuracy_list = []
    if num_runs == 'max':
        features_counts = [i for i in range(1, len(features_rating) + 1)]
    else:
        features_counts = [i for i in range(1, num_runs + 1)]

    for features_count in features_counts:
        if features_counts.index(features_count) % 10 == 0:
            print('Random forest #', str(features_counts.index(features_count)))
        curr_features = features_rating[:features_count]
        curr_features_ids = [positions.index(item) for item in curr_features]
        curr_table = table[:, curr_features_ids].copy()

        clf = RandomForestClassifier(n_estimators=500)
        output = cross_validate(clf, curr_table, y, cv=5, scoring='accuracy', return_estimator=True)
        accuracy_list.append(np.mean(output['test_score']))
    top_accuracy = max(accuracy_list)
    features_top = features_rating[: features_counts[accuracy_list.index(top_accuracy)]]
    return top_accuracy, features_top, accuracy_list, features_rating


def save_results(path, filename, data):
    f = open(path + filename + '.txt', 'w')
    for item in data:
        f.write(str(item) + '\n')
    f.close()


def read_results(path, filename):
    data = []
    f = open(path + filename, 'r')
    for line in f:
        line = line.rstrip()
        data.append(float(line))
    f.close()
    return data
