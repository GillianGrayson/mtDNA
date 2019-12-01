import numpy as np
from tqdm import tqdm
from collections import Counter
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_validate
from bisect import bisect
import itertools


def create_df_freq(raw_data):
    num_nuc = len(raw_data[0][0])
    num_persons = int(np.sum([len(raw_data[i]) for i in range(0, len(raw_data))]))
    df = np.empty(shape=(num_persons, num_nuc), dtype=np.int)
    positions = []

    for nuc_id in tqdm(range(0, num_nuc)):
        curr_nuc = [raw_data[i][j][nuc_id] for i in range(0, len(raw_data)) for j in range(0, len(raw_data[i]))]
        count_dict = Counter(curr_nuc).most_common()
        count_dict = dict(count_dict)
        if '-' in count_dict:
            if count_dict['-'] / len(curr_nuc) > 0.9:
                continue
        reverse_dict = {}
        index_dict = {k: -1 for k in list(count_dict.keys())}
        for key in count_dict:
            if count_dict[key] not in reverse_dict:
                reverse_dict[count_dict[key]] = []
            reverse_dict[count_dict[key]].append(key)
        for key in count_dict:
            if len(reverse_dict[count_dict[key]]) == 1:
                index_dict[key] = list(count_dict.keys()).index(key)
            else:
                for item in reverse_dict[count_dict[key]]:
                    if index_dict[item] == -1:
                        index_dict[item] = min(
                            [list(count_dict.keys()).index(i) for i in reverse_dict[count_dict[key]]])
        for person_id in range(0, num_persons):
            df[person_id, nuc_id] = index_dict[curr_nuc[person_id]]
        if len(np.unique(df[:, nuc_id])) > 1:
            positions.append(nuc_id)

    df = df[:, ~np.all(df[1:] == df[:-1], axis=0)]
    return df, positions


def create_df(raw_data):
    num_nuc = len(raw_data[0][0])
    num_persons = int(np.sum([len(raw_data[i]) for i in range(0, len(raw_data))]))
    df = np.empty(shape=(num_persons, num_nuc), dtype=np.int)
    positions = []

    for nuc_id in tqdm(range(0, num_nuc)):
        curr_nuc = [raw_data[i][j][nuc_id] for i in range(0, len(raw_data)) for j in range(0, len(raw_data[i]))]
        count_dict = Counter(curr_nuc).most_common()
        count_dict = dict(count_dict)
        if '-' in count_dict:
            if count_dict['-'] / len(curr_nuc) > 0.9:
                continue
        if len(count_dict) == 1:
            for person_id in range(0, num_persons):
                df[person_id, nuc_id] = 0
        else:
            for person_id in range(0, num_persons):
                if list(count_dict.keys()).index(curr_nuc[person_id]) == 0:
                    df[person_id, nuc_id] = 0
                else:
                    df[person_id, nuc_id] = 1
        if len(np.unique(df[:, nuc_id])) > 1:
            positions.append(nuc_id)

    df = df[:, ~np.all(df[1:] == df[:-1], axis=0)]
    return df, positions


def run_random_forest(df, classes, positions):
    factor = pd.factorize(classes)
    y = factor[0]
    clf = RandomForestClassifier(n_estimators=500)
    output = cross_validate(clf, df, y, cv=5, scoring='accuracy', return_estimator=True)
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
    return accuracy, features_dict


def run_sequential_random_forest(df, classes, positions):
    factor = pd.factorize(classes)
    y = factor[0]
    clf = RandomForestClassifier(n_estimators=500)
    output = cross_validate(clf, df, y, cv=10, scoring='accuracy', return_estimator=True)

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

    top_features = list(features_dict.keys())
    accuracy_list = []
    features_counts = np.geomspace(2.0, len(top_features) // 10, 10, endpoint=True)
    features_counts = list(set([int(item) for item in features_counts]))
    features_counts.sort()
    for features_count in features_counts:
        if features_counts.index(features_count) % 10 == 0:
            print('Random forest #', str(features_counts.index(features_count)))
        curr_features = top_features[:features_count]
        curr_features_ids = [positions.index(item) for item in curr_features]
        curr_df = df[:, curr_features_ids].copy()

        clf = RandomForestClassifier(n_estimators=500)
        output = cross_validate(clf, curr_df, y, cv=5, scoring='accuracy', return_estimator=True)
        accuracy_list.append(np.mean(output['test_score']))
    top_accuracy = max(accuracy_list)
    top_features = top_features[: features_counts[accuracy_list.index(top_accuracy)]]
    return top_accuracy, top_features


def run_sequential_random_forest_preset(df, classes, positions, features_dict):
    factor = pd.factorize(classes)
    y = factor[0]
    top_features = list(features_dict.keys())
    accuracy_list = []
    features_counts = np.geomspace(2.0, len(top_features), 10, endpoint=True)
    features_counts = list(set([int(item) for item in features_counts]))
    features_counts.sort()
    for features_count in features_counts:
        if features_counts.index(features_count) % 10 == 0:
            print('Random forest #', str(features_counts.index(features_count)))
        curr_features = top_features[:features_count]
        curr_features_ids = [positions.index(item) for item in curr_features]
        curr_df = df[:, curr_features_ids].copy()

        clf = RandomForestClassifier(n_estimators=500)
        output = cross_validate(clf, curr_df, y, cv=5, scoring='accuracy', return_estimator=True)
        accuracy_list.append(np.mean(output['test_score']))
    top_accuracy = max(accuracy_list)
    top_features = top_features[: features_counts[accuracy_list.index(top_accuracy)]]
    return top_accuracy, top_features


def create_features_top(features_dict):
    features = Counter(features_dict)
    top_features = dict([item for item in features.most_common(len(features_dict) // 4)])
    return top_features


def create_regions_stat(top_features, regions):
    top_regions = []
    for feature in top_features:
        position = bisect(regions['start'], feature) - 1
        if feature <= regions['finish'][position]:
            top_regions.append(regions['region'][position])
        else:
            top_regions.append('NA')
    top_regions_dict = dict(Counter(top_regions).most_common())
    regions_sum = sum(top_regions_dict.values(), 0.0)
    for feature in top_regions_dict:
        top_regions_dict[feature] /= regions_sum
    return top_regions_dict


def create_pair_regions_stat(top_features, regions):
    top_regions = []
    for feature_pair in top_features:
        features = [int(feature) for feature in feature_pair.split('_')]
        positions = [bisect(regions['start'], feature) - 1 for feature in features]
        curr_regions = []
        for position in positions:
            if position <= regions['finish'][position]:
                curr_regions.append(regions['region'][position])
            else:
                curr_regions.append('NA')
        curr_regions.sort()
        top_regions.append(' '.join(curr_regions))
    top_regions_dict = dict(Counter(top_regions).most_common())
    regions_sum = sum(top_regions_dict.values(), 0.0)
    features_count = 0
    for feature_pair in list(top_regions_dict.keys()):
        top_regions_dict[feature_pair] /= regions_sum
        if features_count > 40:
            del top_regions_dict[feature_pair]
        features_count += 1
    return top_regions_dict


def get_variable_positions(raw_data):
    num_nuc = len(raw_data[0][0])
    positions = []
    for nuc_id in tqdm(range(0, num_nuc)):
        curr_nuc = [raw_data[i][j][nuc_id] for i in range(0, len(raw_data)) for j in range(0, len(raw_data[i]))]
        count_dict = Counter(curr_nuc).most_common()
        count_dict = dict(count_dict)
        if '-' in count_dict and count_dict['-'] / len(curr_nuc) > 0.9:
            continue
        if len(count_dict) > 1:
            positions.append(nuc_id)
    return positions


def create_co_df_freq(raw_data, variable_positions):
    nuc_combinations = list(itertools.combinations(variable_positions, 2))
    num_persons = int(np.sum([len(raw_data[i]) for i in range(0, len(raw_data))]))
    df = np.empty(shape=(num_persons, len(nuc_combinations)), dtype=np.int)
    positions = []
    combination_id = 0
    for subset in tqdm(nuc_combinations):
        nuc_1_id = subset[0]
        nuc_2_id = subset[1]
        curr_nuc_1 = [raw_data[i][j][nuc_1_id] for i in range(0, len(raw_data)) for j in range(0, len(raw_data[i]))]
        curr_nuc_2 = [raw_data[i][j][nuc_2_id] for i in range(0, len(raw_data)) for j in range(0, len(raw_data[i]))]
        curr_nuc_pair = [curr_nuc_1[i] + curr_nuc_2[i] for i in range(0, num_persons)]
        count_dict = Counter(curr_nuc_pair).most_common()
        count_dict = dict(count_dict)
        if '-' in count_dict:
            if count_dict['-'] / len(curr_nuc_pair) > 0.9:
                continue
        reverse_dict = {}
        index_dict = {k: -1 for k in list(count_dict.keys())}
        for key in count_dict:
            if count_dict[key] not in reverse_dict:
                reverse_dict[count_dict[key]] = []
            reverse_dict[count_dict[key]].append(key)
        for key in count_dict:
            if len(reverse_dict[count_dict[key]]) == 1:
                index_dict[key] = list(count_dict.keys()).index(key)
            else:
                for item in reverse_dict[count_dict[key]]:
                    if index_dict[item] == -1:
                        index_dict[item] = min(
                            [list(count_dict.keys()).index(i) for i in reverse_dict[count_dict[key]]])
        for person_id in range(0, num_persons):
            df[person_id, combination_id] = index_dict[curr_nuc_pair[person_id]]
        if len(np.unique(df[:, combination_id])) > 1:
            positions.append(str(nuc_1_id) + '_' + str(nuc_2_id))
        combination_id += 1

    df = df[:, ~np.all(df[1:] == df[:-1], axis=0)]
    return df, positions


def create_co_df(raw_data, variable_positions):
    nuc_combinations = list(itertools.combinations(variable_positions, 2))
    num_persons = int(np.sum([len(raw_data[i]) for i in range(0, len(raw_data))]))
    df = np.empty(shape=(num_persons, len(nuc_combinations)), dtype=np.int)
    positions = []
    combination_id = 0
    for subset in tqdm(nuc_combinations):
        nuc_1_id = subset[0]
        nuc_2_id = subset[1]
        curr_nuc_1 = [raw_data[i][j][nuc_1_id] for i in range(0, len(raw_data)) for j in range(0, len(raw_data[i]))]
        curr_nuc_2 = [raw_data[i][j][nuc_2_id] for i in range(0, len(raw_data)) for j in range(0, len(raw_data[i]))]
        curr_nuc_pair = [curr_nuc_1[i] + curr_nuc_2[i] for i in range(0, num_persons)]
        count_dict = Counter(curr_nuc_pair).most_common()
        count_dict = dict(count_dict)
        if '-' in count_dict:
            if count_dict['-'] / len(curr_nuc_pair) > 0.9:
                continue
        if len(count_dict) == 1:
            for person_id in range(0, num_persons):
                df[person_id, combination_id] = 0
        else:
            for person_id in range(0, num_persons):
                if list(count_dict.keys()).index(curr_nuc_pair[person_id]) == 0:
                    df[person_id, combination_id] = 0
                else:
                    df[person_id, combination_id] = 1

        if len(np.unique(df[:, combination_id])) > 1:
            positions.append(str(nuc_1_id) + '_' + str(nuc_2_id))

        combination_id += 1

    df = df[:, ~np.all(df[1:] == df[:-1], axis=0)]
    return df, positions


def get_haplogroups_statistics(haplogroups):
    haplogroups_stat = {}
    for group in haplogroups:
        haplogroups_stat[group] = dict(Counter(haplogroups[group]).most_common())
    return haplogroups_stat
