import numpy as np
import pandas as pd
from tqdm import tqdm
from collections import Counter
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_validate
from bisect import bisect


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


def read_haplogroups(data_path, classes):
    data = pd.read_excel(data_path + 'subjects.xlsx').to_dict('list')
    data_classes = [item for curr_class in classes for item in classes[curr_class]]
    haplogroups = []
    for item_id in range(0, len(data['subject'])):
        if data['height'][item_id] in data_classes:
            haplogroup = data['group'][item_id]
            if haplogroup not in haplogroups:
                haplogroups.append(haplogroup)
    return haplogroups


def get_haplogroups_positions(data_path, haplogroups):
    phylotrees = pd.read_excel(data_path + 'phylotrees.xlsx').to_dict('list')
    positions = []
    for item_id in tqdm(range(0, len(phylotrees['phylotree']))):
        if phylotrees['haplogroup'][item_id] in haplogroups:
            position = phylotrees['position'][item_id]
            if isinstance(position, str):
                start_position = int(position.split('-')[0])
                end_position = int(position.split('-')[1])
                curr_positions = list(range(start_position, end_position + 1))
                for curr_position in curr_positions:
                    if curr_position not in positions:
                        positions.append(curr_position)
            else:
                if position not in positions:
                    positions.append(position)
    positions.sort()
    return positions


def remove_items_from_list(initial_list, positions_to_remove):
    modified_list = initial_list.copy()
    for item in positions_to_remove:
        if item in initial_list:
            modified_list.remove(item)
    return modified_list


def calculate_mutation_frequency(data, classes, features):
    frequency_dict = {str(position): {data_class: 0.0 for data_class in classes} for position in features}
    for position in features:
        for data_class in classes:
            data_class_id = classes.index(data_class)
            curr_nuc = [data[data_class_id][person_id][position] for person_id in range(0, len(data[data_class_id]))]
            count_dict = Counter(curr_nuc).most_common()
            count_dict = dict(count_dict)
            count_list = [count_dict[item] for item in count_dict]
            frequency_dict[str(position)][data_class] = np.sum(count_list[1:]) / np.sum(count_list)
    return frequency_dict


def filter_frequency_dict(frequency_dict):
    positions = list(frequency_dict.keys())
    frequency_dict_1 = {}
    frequency_dict_2 = {}
    frequency_dict_3 = {}
    for position in positions:
        populations = list(frequency_dict[position].keys())
        is_all_non_zero = 0
        for population in populations:
            if frequency_dict[position][population] > 0:
                is_all_non_zero += 1
        if is_all_non_zero == 1:
            frequency_dict_1[position] = {population: frequency_dict[position][population] for population in
                                          populations}
        if is_all_non_zero == 2:
            frequency_dict_2[position] = {population: frequency_dict[position][population] for population in
                                          populations}
        if is_all_non_zero == 3:
            frequency_dict_3[position] = {population: frequency_dict[position][population] for population in
                                          populations}
    return frequency_dict_1, frequency_dict_2, frequency_dict_3


def calculate_regions_statistics(features, regions):
    regions_stat = []
    for feature in features:
        if isinstance(feature, str):
            feature = int(feature)
        position = bisect(regions['start'], feature) - 1
        if feature <= regions['finish'][position]:
            regions_stat.append(regions['region'][position])
        else:
            regions_stat.append('NA')
    regions_dict = dict(Counter(regions_stat).most_common())
    regions_sum = sum(regions_dict.values(), 0.0)
    for feature in regions_dict:
        regions_dict[feature] /= regions_sum
    return regions_dict


def create_mutation_statistics(data, positions, classes):
    result_dict = {'Position': []}
    classes_extended = []
    for curr_class in classes:
        classes_extended.append(curr_class + ' Main')
        classes_extended.append(curr_class + ' Minor')
        classes_extended.append(curr_class + ' Other')
    result_dict.update({curr_class: [] for curr_class in classes_extended})
    for position in positions:
        result_dict['Position'].append(position + 1)
        for curr_class in classes:
            curr_nucleotide = [data[curr_class][person_id][position] for person_id in range(0, len(data[curr_class]))]
            count_dict = Counter(curr_nucleotide).most_common()
            count_dict = dict(count_dict)
            variants_list = list(count_dict.keys())
            result_dict[curr_class + ' Main'].append(variants_list[0])
            if len(variants_list) == 1:
                result_dict[curr_class + ' Minor'].append('')
                result_dict[curr_class + ' Other'].append('')
            if len(variants_list) > 1:
                result_dict[curr_class + ' Minor'].append(variants_list[1])
                if len(variants_list) == 2:
                    result_dict[curr_class + ' Other'].append(''.join(variants_list[2:]))
            if len(variants_list) > 2:
                result_dict[curr_class + ' Other'].append(''.join(variants_list[2:]))
    return result_dict
