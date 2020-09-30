import numpy as np
import pandas as pd
from tqdm import tqdm
from collections import Counter
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_validate
from bisect import bisect
import itertools


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


def get_variable_positions(raw_data):
    data_classes = list(raw_data.keys())
    num_nucleotides = len(raw_data[data_classes[0]][0])
    variable_positions = []
    for nuc_id in tqdm(range(0, num_nucleotides)):
        curr_nuc = [raw_data[data_class][j][nuc_id] for data_class in data_classes for j in
                    range(0, len(raw_data[data_class]))]
        count_dict = Counter(curr_nuc).most_common()
        count_dict = dict(count_dict)
        if '-' in count_dict and count_dict['-'] / len(curr_nuc) > 0.9:
            continue
        if len(count_dict) > 1:
            variable_positions.append(nuc_id)
    return variable_positions


def create_classes_table_pairs(raw_data, variable_positions):
    data_classes = list(raw_data.keys())
    nucleotide_pairs = list(itertools.combinations(variable_positions, 2))
    num_persons = int(np.sum([len(raw_data[data_class]) for data_class in data_classes]))
    table = np.empty(shape=(num_persons, len(nucleotide_pairs)), dtype=np.int)
    mut_positions = []
    combination_id = 0
    for nucleotide_pair in tqdm(nucleotide_pairs):
        nuc_1_id = nucleotide_pair[0]
        nuc_2_id = nucleotide_pair[1]
        curr_nuc_1 = [raw_data[data_class][person_id][nuc_1_id] for data_class in data_classes for person_id in
                      range(0, len(raw_data[data_class]))]
        curr_nuc_2 = [raw_data[data_class][person_id][nuc_2_id] for data_class in data_classes for person_id in
                      range(0, len(raw_data[data_class]))]
        curr_nuc_pair = [curr_nuc_1[i] + curr_nuc_2[i] for i in range(0, num_persons)]
        count_dict = Counter(curr_nuc_pair).most_common()
        count_dict = dict(count_dict)
        if '--' in count_dict:
            if count_dict['--'] / len(curr_nuc_pair) > 0.9:
                continue
        if len(count_dict) == 1:
            continue
        else:
            for person_id in range(0, num_persons):
                if list(count_dict.keys()).index(curr_nuc_pair[person_id]) == 0:
                    table[person_id, combination_id] = 0
                else:
                    table[person_id, combination_id] = 1
        if len(np.unique(table[:, combination_id])) > 1:
            mut_positions.append(nucleotide_pair)
        combination_id += 1
    table = table[:, ~np.all(table[1:] == table[:-1], axis=0)]
    return table, mut_positions


def run_sequential_random_forest(table, classes, positions, num_runs):
    factor = pd.factorize(classes)
    y = factor[0]
    clf = RandomForestClassifier(n_estimators=500)
    output = cross_validate(clf, table, y, cv=10, scoring='accuracy', return_estimator=True)
    accuracy = np.mean(output['test_score'])

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
    accuracy_list = [accuracy]
    if num_runs == 'max':
        features_counts = [i for i in range(1, len(features_rating) + 1)]
    elif isinstance(num_runs, str) and num_runs.startswith('>'):
        features_counts = [1]
    else:
        features_counts = [i for i in range(1, num_runs + 1)]

    for features_count in features_counts:
        if features_counts.index(features_count) % 10 == 0:
            print('Random forest #', str(features_counts.index(features_count)))
        curr_features = features_rating[:features_count]
        curr_features_ids = [positions.index(item) for item in curr_features]
        curr_table = table[:, curr_features_ids].copy()

        clf = RandomForestClassifier(n_estimators=500)
        output = cross_validate(clf, curr_table, y, cv=10, scoring='accuracy', return_estimator=True)
        accuracy_list.append(np.mean(output['test_score']))
    top_accuracy = max(accuracy_list)
    top_index = accuracy_list.index(top_accuracy)
    if top_index == 0:
        if isinstance(num_runs, str) and num_runs.startswith('>'):
            target_value = float(num_runs[1:])
            features_top = []
            for feature in features_rating:
                if features_dict[feature] > target_value:
                    features_top.append(feature)
        else:
            features_top = features_rating
    else:
        features_top = features_rating[: top_index]
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


def remove_items_from_pair(initial_list, positions_to_remove):
    modified_list = initial_list.copy()
    for item in initial_list:
        if item[0] in positions_to_remove or item[1] in positions_to_remove:
            modified_list.remove(item)
    return modified_list


def remove_pairs_from_list(initial_list, positions_to_remove):
    modified_list = initial_list.copy()
    for item in initial_list:
        if item[0] in positions_to_remove and item[1] in positions_to_remove:
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
        classes_extended.append(curr_class + ' Main Freq')
        classes_extended.append(curr_class + ' Minor')
        classes_extended.append(curr_class + ' Minor Freq')
        classes_extended.append(curr_class + ' Other')
        classes_extended.append(curr_class + ' Other Freq')
    result_dict.update({curr_class: [] for curr_class in classes_extended})
    for position in positions:
        result_dict['Position'].append(position + 1)
        for curr_class in classes:
            curr_nucleotide = [data[curr_class][person_id][position] for person_id in range(0, len(data[curr_class]))]
            count_dict = Counter(curr_nucleotide).most_common()
            count_dict = dict(count_dict)
            variants_list = list(count_dict.keys())
            result_dict[curr_class + ' Main'].append(variants_list[0])
            result_dict[curr_class + ' Main Freq'].append(count_dict[variants_list[0]] / len(curr_nucleotide))
            if len(variants_list) == 1:
                result_dict[curr_class + ' Minor'].append('')
                result_dict[curr_class + ' Minor Freq'].append(0.0)
                result_dict[curr_class + ' Other'].append('')
                result_dict[curr_class + ' Other Freq'].append(0.0)
            if len(variants_list) > 1:
                result_dict[curr_class + ' Minor'].append(variants_list[1])
                result_dict[curr_class + ' Minor Freq'].append(count_dict[variants_list[1]] / len(curr_nucleotide))
                if len(variants_list) == 2:
                    result_dict[curr_class + ' Other'].append('')
                    result_dict[curr_class + ' Other Freq'].append(0.0)
            if len(variants_list) > 2:
                result_dict[curr_class + ' Other'].append(''.join(variants_list[2:]))
                result_dict[curr_class + ' Other Freq'].append(
                    np.sum([count_dict[variants_list[i]] for i in range(2, len(variants_list))]) / len(curr_nucleotide))
    return result_dict


def create_pair_statistics(data, positions, classes):
    result_dict = {'Position': []}
    classes_extended = []
    for curr_class in classes:
        classes_extended.append(curr_class + ' Main')
        classes_extended.append(curr_class + ' Main Freq')
        classes_extended.append(curr_class + ' Minor')
        classes_extended.append(curr_class + ' Minor Freq')
        classes_extended.append(curr_class + ' Other')
        classes_extended.append(curr_class + ' Other Freq')
    result_dict.update({curr_class: [] for curr_class in classes_extended})
    for position_pair in tqdm(positions):
        if isinstance(position_pair, str):
            position_1 = int(position_pair[1:-1].split(',')[0])
            position_2 = int(position_pair[1:-1].split(',')[1])
        else:
            position_1 = int(position_pair[0])
            position_2 = int(position_pair[1])
        result_dict['Position'].append(str(position_1 + 1) + ', ' + str(position_2 + 1))
        for curr_class in classes:
            curr_nucleotide_1 = [data[curr_class][person_id][position_1] for person_id in
                                 range(0, len(data[curr_class]))]
            curr_nucleotide_2 = [data[curr_class][person_id][position_2] for person_id in
                                 range(0, len(data[curr_class]))]
            curr_nuc_pair = [curr_nucleotide_1[i] + curr_nucleotide_2[i] for i in range(0, len(curr_nucleotide_1))]
            count_dict = Counter(curr_nuc_pair).most_common()
            count_dict = dict(count_dict)
            variants_list = list(count_dict.keys())
            result_dict[curr_class + ' Main'].append(variants_list[0])
            result_dict[curr_class + ' Main Freq'].append(count_dict[variants_list[0]] / len(curr_nuc_pair))
            if len(variants_list) == 1:
                result_dict[curr_class + ' Minor'].append('')
                result_dict[curr_class + ' Minor Freq'].append(0.0)
                result_dict[curr_class + ' Other'].append('')
                result_dict[curr_class + ' Other Freq'].append(0.0)
            if len(variants_list) > 1:
                result_dict[curr_class + ' Minor'].append(variants_list[1])
                result_dict[curr_class + ' Minor Freq'].append(count_dict[variants_list[1]] / len(curr_nuc_pair))
                if len(variants_list) == 2:
                    result_dict[curr_class + ' Other'].append('')
                    result_dict[curr_class + ' Other Freq'].append(0.0)
            if len(variants_list) > 2:
                result_dict[curr_class + ' Other'].append(''.join(variants_list[2:]))
                result_dict[curr_class + ' Other Freq'].append(
                    np.sum([count_dict[variants_list[i]] for i in range(2, len(variants_list))]) / len(curr_nuc_pair))
    return result_dict


def create_selected_statistics(data, positions, classes):
    result_dict = {'Position': []}
    classes_extended = []
    for curr_class in classes:
        classes_extended.append(curr_class + ' Main')
        classes_extended.append(curr_class + ' Main Freq')
        classes_extended.append(curr_class + ' Minor')
        classes_extended.append(curr_class + ' Minor Freq')
        classes_extended.append(curr_class + ' Other')
        classes_extended.append(curr_class + ' Other Freq')
    result_dict.update({curr_class: [] for curr_class in classes_extended})
    for position_pair in tqdm(positions):
        curr_dict = {'Position': ''}
        curr_dict.update({curr_class: 0.0 for curr_class in classes_extended})
        if isinstance(position_pair, str):
            position_1 = int(position_pair[1:-1].split(',')[0])
            position_2 = int(position_pair[1:-1].split(',')[1])
        else:
            position_1 = int(position_pair[0])
            position_2 = int(position_pair[1])
        curr_dict['Position'] = str(position_1 + 1) + ', ' + str(position_2 + 1)
        for curr_class in classes:
            curr_nucleotide_1 = [data[curr_class][person_id][position_1] for person_id in
                                 range(0, len(data[curr_class]))]
            curr_nucleotide_2 = [data[curr_class][person_id][position_2] for person_id in
                                 range(0, len(data[curr_class]))]
            curr_nuc_pair = [curr_nucleotide_1[i] + curr_nucleotide_2[i] for i in range(0, len(curr_nucleotide_1))]
            count_dict = Counter(curr_nuc_pair).most_common()
            count_dict = dict(count_dict)
            variants_list = list(count_dict.keys())
            curr_dict[curr_class + ' Main'] = variants_list[0]
            curr_dict[curr_class + ' Main Freq'] = count_dict[variants_list[0]] / len(curr_nuc_pair)
            if len(variants_list) == 1:
                curr_dict[curr_class + ' Minor'] = ''
                curr_dict[curr_class + ' Minor Freq'] = 0.0
                curr_dict[curr_class + ' Other'] = ''
                curr_dict[curr_class + ' Other Freq'] = 0.0
            if len(variants_list) > 1:
                curr_dict[curr_class + ' Minor'] = variants_list[1]
                curr_dict[curr_class + ' Minor Freq'] = count_dict[variants_list[1]] / len(curr_nuc_pair)
                if len(variants_list) == 2:
                    curr_dict[curr_class + ' Other'] = ''
                    curr_dict[curr_class + ' Other Freq'] = 0.0
            if len(variants_list) > 2:
                curr_dict[curr_class + ' Other'] = ''.join(variants_list[2:])
                curr_dict[curr_class + ' Other Freq'] = np.sum(
                    [count_dict[variants_list[i]] for i in range(2, len(variants_list))]) / len(curr_nuc_pair)

        if curr_dict['Asian Low Altitude Main Freq'] > curr_dict['Tibetan High Altitude Main Freq'] and \
                curr_dict['Asian Low Altitude Main Freq'] > curr_dict['Tibetan Main Freq'] and \
                curr_dict['Asian Low Altitude Main Freq'] > curr_dict['Andes Main Freq'] and \
                curr_dict['Asian Low Altitude Main Freq'] > curr_dict['Ethiopia Main Freq']:
            result_dict['Position'].append(curr_dict['Position'])
            for curr_class in classes:
                result_dict[curr_class + ' Main'].append(curr_dict[curr_class + ' Main'])
                result_dict[curr_class + ' Main Freq'].append(curr_dict[curr_class + ' Main Freq'])
                result_dict[curr_class + ' Minor'].append(curr_dict[curr_class + ' Minor'])
                result_dict[curr_class + ' Minor Freq'].append(curr_dict[curr_class + ' Minor Freq'])
                result_dict[curr_class + ' Other'].append(curr_dict[curr_class + ' Other'])
                result_dict[curr_class + ' Other Freq'].append(curr_dict[curr_class + ' Other Freq'])
        elif curr_dict['Asian Low Altitude Main Freq'] < curr_dict['Tibetan High Altitude Main Freq'] and \
                curr_dict['Asian Low Altitude Main Freq'] < curr_dict['Tibetan Main Freq'] and \
                curr_dict['Asian Low Altitude Main Freq'] < curr_dict['Andes Main Freq'] and \
                curr_dict['Asian Low Altitude Main Freq'] < curr_dict['Ethiopia Main Freq']:
            result_dict['Position'].append(curr_dict['Position'])
            for curr_class in classes:
                result_dict[curr_class + ' Main'].append(curr_dict[curr_class + ' Main'])
                result_dict[curr_class + ' Main Freq'].append(curr_dict[curr_class + ' Main Freq'])
                result_dict[curr_class + ' Minor'].append(curr_dict[curr_class + ' Minor'])
                result_dict[curr_class + ' Minor Freq'].append(curr_dict[curr_class + ' Minor Freq'])
                result_dict[curr_class + ' Other'].append(curr_dict[curr_class + ' Other'])
                result_dict[curr_class + ' Other Freq'].append(curr_dict[curr_class + ' Other Freq'])
    return result_dict


def calculate_pair_regions_statistics(features, regions):
    regions_stat = []
    for position_pair in tqdm(features):
        if isinstance(position_pair, str):
            feature_1 = int(position_pair[1:-1].split(',')[0])
            feature_2 = int(position_pair[1:-1].split(',')[1])
        else:
            feature_1 = int(position_pair[0])
            feature_2 = int(position_pair[1])
        position_1 = bisect(regions['start'], feature_1) - 1
        if feature_1 <= regions['finish'][position_1]:
            region_1 = regions['region'][position_1]
        else:
            region_1 = 'NA'
        position_2 = bisect(regions['start'], feature_2) - 1
        if feature_2 <= regions['finish'][position_2]:
            region_2 = regions['region'][position_2]
        else:
            region_2 = 'NA'
        curr_regions = [region_1, region_2]
        curr_regions.sort()
        regions_stat.append(';'.join(curr_regions))
    regions_dict = dict(Counter(regions_stat).most_common())
    regions_sum = sum(regions_dict.values(), 0.0)
    for feature in regions_dict:
        regions_dict[feature] /= regions_sum
    result_dict = {'Region': [], 'Freq': []}
    for key in regions_dict:
        result_dict['Region'].append(key)
        result_dict['Freq'].append(regions_dict[key])
    return result_dict
