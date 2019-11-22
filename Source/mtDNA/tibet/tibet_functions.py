import numpy as np
from tqdm import tqdm
from collections import Counter
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_validate
from bisect import bisect
import plotly
import plotly.graph_objs as go


def create_df_freq(raw_data):

    num_nuc = len(raw_data[0][0])
    num_persons = int(np.sum([len(raw_data[i]) for i in range(0, len(raw_data))]))
    df = np.empty(shape=(num_persons, num_nuc), dtype=np.int)
    positions = []

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
    return df


def run_random_forest(df, classes, positions):
    factor = pd.factorize(classes)
    y = factor[0]
    clf = RandomForestClassifier(n_estimators=500)
    output = cross_validate(clf, df, y, cv=10, scoring='accuracy', return_estimator=True)
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


def create_features_top(features_dict):
    features = Counter(features_dict)
    top_features = dict([item for item in features.most_common(len(features_dict)//4)])
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


def plot_hist(data, suffix, file_path):
    fig = go.Figure(go.Bar(
        x=list(data.values())[::-1],
        y=list(data.keys())[::-1],
        orientation='h'
    ))
    fig.update_yaxes(
        tickfont=dict(size=10)
    )
    fig.update_layout(width=700,
                      height=1000)

    plotly.offline.plot(fig, filename=file_path + suffix + '_hist.html', auto_open=False, show_link=True)
    plotly.io.write_image(fig, file_path + suffix + '_hist.png')
    plotly.io.write_image(fig, file_path + suffix + '_hist.pdf')


def save_dict(data, filename):
    f = open(filename, 'w')
    for key in data.keys():
        f.write(str(key) + '\t' + str(data[key]) + '\n')
    f.close()


def save_list(data, filename):
    f = open(filename, 'w')
    for item in data:
        f.write(str(item) + '\n')
    f.close()
