import numpy as np
from tqdm import tqdm
from collections import Counter
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_validate


def create_df_freq(raw_data):

    num_nuc = len(raw_data[0][0])
    num_persons = int(np.sum([len(raw_data[i]) for i in range(0, len(raw_data))]))
    df = np.empty(shape=(num_persons, num_nuc), dtype=np.int)

    for nuc_id in tqdm(range(0, num_nuc)):
        curr_nuc = [raw_data[i][j][nuc_id] for i in range(0, len(raw_data)) for j in range(0, len(raw_data[i]))]
        count_dict = Counter(curr_nuc).most_common()
        count_dict = dict(count_dict)
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

    df = df[:, ~np.all(df[1:] == df[:-1], axis=0)]
    return df


def create_df(raw_data):

    num_nuc = len(raw_data[0][0])
    num_persons = int(np.sum([len(raw_data[i]) for i in range(0, len(raw_data))]))
    df = np.empty(shape=(num_persons, num_nuc), dtype=np.int)

    for nuc_id in tqdm(range(0, num_nuc)):
        curr_nuc = [raw_data[i][j][nuc_id] for i in range(0, len(raw_data)) for j in range(0, len(raw_data[i]))]
        count_dict = Counter(curr_nuc).most_common()
        count_dict = dict(count_dict)
        if len(count_dict) == 1:
            for person_id in range(0, num_persons):
                df[person_id, nuc_id] = 0
        else:
            for person_id in range(0, num_persons):
                if list(count_dict.keys()).index(curr_nuc[person_id]) == 0:
                    df[person_id, nuc_id] = 0
                else:
                    df[person_id, nuc_id] = 1

    df = df[:, ~np.all(df[1:] == df[:-1], axis=0)]
    return df


def run_random_forest(df, classes):
    factor = pd.factorize(classes)
    y = factor[0]
    clf = RandomForestClassifier(n_estimators=500)
    output = cross_validate(clf, df, y, cv=10, scoring='accuracy', return_estimator=True)
    accuracy = np.mean(output['test_score'])
    return accuracy
