from save import save_df
import pandas as pd
import numpy as np
import random
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_validate
import time


# type 0 - calculating of 2 reference frequencies
def rf_type_0_mt(config, results):
    reference_part = float(config.params_dict['reference_part'])
    reference_pop = config.params_dict['reference_pop']
    target_pop = config.params_dict['target_pop']
    reference_size = int(len(config.pop_person_dict[reference_pop]) * reference_part)
    reference_list = random.sample(config.pop_person_dict[reference_pop], reference_size)
    reference_frequencies = [0, 0]

    for task in config.params_dict['config_mt_genes']:

        genes_ids = list(task)
        genes_names = [config.params_dict['genes_list'][0][i] for i in genes_ids]

        print(';'.join(genes_names))

        # Reference group

        persons = config.person_index_dict

        target_samples_ids_mt = []
        target_samples_names_mt = []

        for sample_name in persons:
            if sample_name in reference_list:
                target_samples_ids_mt.append(persons[sample_name])
                target_samples_names_mt.append(sample_name)

        number_mt_snps = 0

        if int(config.params_dict['run_timer']) == 1:
            start_ref = time.process_time()

        for gene_id in range(0, len(genes_ids)):
            gene_index = config.data_position_dict[genes_names[gene_id]]
            for row in config.data[gene_index]:
                snp_data_mt = list(row[i] for i in target_samples_ids_mt)
                for i in range(0, len(snp_data_mt)):
                    if snp_data_mt[i] == 0:
                        reference_frequencies[0] += 1
                    elif snp_data_mt[i] == 1:
                        reference_frequencies[1] += 1
                number_mt_snps += 1

        reference_frequencies = [freq / sum(reference_frequencies) for freq in reference_frequencies]

        if int(config.params_dict['run_timer']) == 1:
            print('Time for frequencies calculating: ' + str(time.process_time() - start_ref))

        # Remaining group

        target_samples_ids_mt = []
        target_samples_names_mt = []

        for sample_name in persons:
            if sample_name in config.pop_person_dict[target_pop] \
                    or sample_name in config.pop_person_dict[reference_pop]:
                target_samples_ids_mt.append(persons[sample_name])
                target_samples_names_mt.append(sample_name)

        df_ref_mt = np.empty(shape=(len(target_samples_names_mt), number_mt_snps), dtype=np.float32)
        names_mt = []

        line_count_mt = 0

        if int(config.params_dict['run_timer']) == 1:
            start_df = time.process_time()

        for gene_id in range(0, len(genes_ids)):
            gene_index = config.data_position_dict[genes_names[gene_id]]
            row_id = 0
            for row in config.data[gene_index]:
                snp_data_mt = list(row[i] for i in target_samples_ids_mt)
                if len(set(snp_data_mt)) == 1:
                    continue

                combination_data = []

                for i in range(0, len(snp_data_mt)):
                    if snp_data_mt[i] == 0:
                        combination_data.append(1 - reference_frequencies[0])
                    elif snp_data_mt[i] == 1:
                        combination_data.append(1 - reference_frequencies[1])

                df_ref_mt[:, line_count_mt] = combination_data

                if len(set(combination_data)) > 1:
                    gene_mt = config.params_dict['genes_list'][0][genes_ids[gene_id]]
                    snp_pos_mt = [name for name, index in config.gene_snp_dict[gene_mt].items() if index == row_id][0]
                    name_mt = gene_mt + '_' + snp_pos_mt
                    if name_mt not in names_mt:
                        names_mt.append(name_mt)

                line_count_mt += 1
                row_id += 1

        if int(config.params_dict['run_timer']) == 1:
            print('Time for data frame creating: ' + str(time.process_time() - start_df))

        df_ref_mt = df_ref_mt[:, : len(names_mt)]

        data_classes = []
        for item in target_samples_names_mt:
            if item in config.pop_person_dict[target_pop]:
                data_classes.append(target_pop)
            elif item in config.pop_person_dict[reference_pop]:
                data_classes.append(reference_pop)

        factor = pd.factorize(data_classes)
        y = factor[0]

        if int(config.params_dict['run_timer']) == 1:
            start_rf = time.process_time()

        clf = RandomForestClassifier(n_estimators=100)
        output = cross_validate(clf, df_ref_mt, y, cv=10, scoring='accuracy', return_estimator=True)
        accuracy = np.mean(output['test_score'])
        if int(config.params_dict['run_timer']) == 1:
            score_time = output['score_time']
            print('Cross validation time: ' + ', '.join([str(item) for item in score_time]))
            print('Total time for random forest ' + str(time.process_time() - start_rf))
            print('Mean accuracy: ' + str(accuracy))

        features_dict = dict((key, []) for key in names_mt)

        num_features = 0

        for idx, estimator in enumerate(output['estimator']):
            feature_importances = pd.DataFrame(estimator.feature_importances_,
                                               index=names_mt,
                                               columns=['importance']).sort_values('importance', ascending=False)

            features_names = list(feature_importances.index.values)
            features_values = list(feature_importances.values)
            for i in range(0, len(features_names)):
                features_dict[features_names[i]].append(features_values[i][0])

        for key in features_dict.keys():
            features_dict[key] = np.mean(features_dict[key])
            if features_dict[key] > 0.001:
                num_features += 1

        if accuracy >= float(config.params_dict['target_accuracy']):
            results.accuracy.append(accuracy)
            results.num_features.append(num_features)
            results.mt_genes.append(genes_ids)

        if int(config.params_dict['num_features']) > 0:
            features_dict = {k: v for k, v in sorted(features_dict.items(), reverse=True, key=lambda x: x[1])}
            results.features = features_dict


# type 1 - using genes variations for rf
def rf_type_1_mt(config, results):
    reference_pop = config.params_dict['reference_pop']
    target_pop = config.params_dict['target_pop']

    for task in config.params_dict['config_mt_genes']:

        genes_ids = list(task)
        genes_names = [config.params_dict['genes_list'][0][i] for i in genes_ids]

        print(';'.join(genes_names))

        persons = config.person_index_dict

        target_samples_ids_mt = []
        target_samples_names_mt = []

        for sample_name in persons:
            if sample_name in config.pop_person_dict[target_pop] \
                    or sample_name in config.pop_person_dict[reference_pop]:
                target_samples_ids_mt.append(persons[sample_name])
                target_samples_names_mt.append(sample_name)

        number_mt_snps = 0
        for gene_id in range(0, len(genes_ids)):
            gene_index = config.data_position_dict[genes_names[gene_id]]
            number_mt_snps += config.data[gene_index].shape[0]

        df_ref_mt = np.empty(shape=(len(target_samples_names_mt), number_mt_snps), dtype=np.int)
        names_mt = []

        line_count_mt = 0

        if int(config.params_dict['run_timer']) == 1:
            start_df = time.process_time()

        for gene_id in range(0, len(genes_ids)):
            gene_index = config.data_position_dict[genes_names[gene_id]]
            row_id = 0
            for row in config.data[gene_index]:
                snp_data_mt = list(row[i] for i in target_samples_ids_mt)
                if len(set(snp_data_mt)) == 1:
                    continue

                combination_data = []

                for i in range(0, len(snp_data_mt)):
                    if snp_data_mt[i] == 0:
                        combination_data.append(0)
                    elif snp_data_mt[i] == 1:
                        combination_data.append(1)

                df_ref_mt[:, line_count_mt] = combination_data

                if len(set(combination_data)) > 1:
                    gene_mt = config.params_dict['genes_list'][0][genes_ids[gene_id]]
                    snp_pos_mt = [name for name, index in config.gene_snp_dict[gene_mt].items() if index == row_id][0]
                    name_mt = gene_mt + '_' + snp_pos_mt
                    if name_mt not in names_mt:
                        names_mt.append(name_mt)

                line_count_mt += 1
                row_id += 1

        if int(config.params_dict['run_timer']) == 1:
            print('Time for data frame creating: ' + str(time.process_time() - start_df))

        df_ref_mt = df_ref_mt[:, : len(names_mt)]

        data_classes = []
        for item in target_samples_names_mt:
            if item in config.pop_person_dict[target_pop]:
                data_classes.append(target_pop)
            elif item in config.pop_person_dict[reference_pop]:
                data_classes.append(reference_pop)

        factor = pd.factorize(data_classes)
        y = factor[0]

        if int(config.params_dict['run_timer']) == 1:
            start_rf = time.process_time()

        clf = RandomForestClassifier(n_estimators=100)
        output = cross_validate(clf, df_ref_mt, y, cv=10, scoring='accuracy', return_estimator=True)
        accuracy = np.mean(output['test_score'])
        if int(config.params_dict['run_timer']) == 1:
            score_time = output['score_time']
            print('Cross validation time: ' + ', '.join([str(item) for item in score_time]))
            print('Total time for random forest ' + str(time.process_time() - start_rf))
            print('Mean accuracy: ' + str(accuracy))

        features_dict = dict((key, []) for key in names_mt)

        num_features = 0

        for idx, estimator in enumerate(output['estimator']):
            feature_importances = pd.DataFrame(estimator.feature_importances_,
                                               index=names_mt,
                                               columns=['importance']).sort_values('importance', ascending=False)

            features_names = list(feature_importances.index.values)
            features_values = list(feature_importances.values)
            for i in range(0, len(features_names)):
                features_dict[features_names[i]].append(features_values[i][0])

        for key in features_dict.keys():
            features_dict[key] = np.mean(features_dict[key])
            if features_dict[key] > 0.001:
                num_features += 1

        if accuracy >= float(config.params_dict['target_accuracy']):
            results.accuracy.append(accuracy)
            results.num_features.append(num_features)
            results.mt_genes.append(genes_ids)

        if int(config.params_dict['num_features']) > 0:
            features_dict = {k: v for k, v in sorted(features_dict.items(), reverse=True, key=lambda x: x[1])}
            results.features = features_dict


# type 2 - using mean frequencies over genes for rf
def rf_type_2_mt(config, results):
    reference_part = float(config.params_dict['reference_part'])
    reference_pop = config.params_dict['reference_pop']
    target_pop = config.params_dict['target_pop']
    reference_size = int(len(config.pop_person_dict[reference_pop]) * reference_part)
    reference_list = random.sample(config.pop_person_dict[reference_pop], reference_size)

    for task in config.params_dict['config_mt_genes']:

        genes_ids = list(task)
        genes_names = [config.params_dict['genes_list'][0][i] for i in genes_ids]

        print(';'.join(genes_names))

        # Reference group

        reference_frequencies = np.zeros((len(genes_ids), 2), dtype=np.float32)

        persons = config.person_index_dict

        target_samples_ids_mt = []
        target_samples_names_mt = []

        for sample_name in persons:
            if sample_name in reference_list:
                target_samples_ids_mt.append(persons[sample_name])
                target_samples_names_mt.append(sample_name)

        if int(config.params_dict['run_timer']) == 1:
            start_ref = time.process_time()

        for gene_id in range(0, len(genes_ids)):
            gene_index = config.data_position_dict[genes_names[gene_id]]
            for row in config.data[gene_index]:
                snp_data_mt = list(row[i] for i in target_samples_ids_mt)
                for i in range(0, len(snp_data_mt)):
                    if snp_data_mt[i] == 0:
                        reference_frequencies[gene_id, 0] += 1
                    elif snp_data_mt[i] == 1:
                        reference_frequencies[gene_id, 1] += 1

        for gene_id in range(0, len(genes_ids)):
            ref_sum = np.sum(reference_frequencies[gene_id, :])
            for i in range(0, 2):
                reference_frequencies[gene_id, i] = reference_frequencies[gene_id, i] / ref_sum

        if int(config.params_dict['run_timer']) == 1:
            print('Time for frequencies calculating: ' + str(time.process_time() - start_ref))

        # Remaining group

        target_samples_ids_mt = []
        target_samples_names_mt = []

        for sample_name in persons:
            if sample_name in config.pop_person_dict[target_pop] \
                    or sample_name in config.pop_person_dict[reference_pop]:
                target_samples_ids_mt.append(persons[sample_name])
                target_samples_names_mt.append(sample_name)

        num_frequencies = len(genes_ids)

        df_ref_mt = np.empty(shape=(len(target_samples_names_mt), num_frequencies), dtype=np.float32)

        if int(config.params_dict['run_timer']) == 1:
            start_df = time.process_time()

        names_mt = []

        for gene_id in range(0, len(genes_ids)):
            gene_index = config.data_position_dict[genes_names[gene_id]]
            gene_data = np.zeros(shape=len(target_samples_ids_mt), dtype=np.float32)
            num_snps = 0
            for row in config.data[gene_index]:
                snp_data_mt = list(row[i] for i in target_samples_ids_mt)
                if len(set(snp_data_mt)) == 1:
                    continue

                for i in range(0, len(snp_data_mt)):
                    if snp_data_mt[i] == 0:
                        gene_data[i] += 1 - reference_frequencies[gene_id, 0]
                    elif snp_data_mt[i] == 1:
                        gene_data[i] += 1 - reference_frequencies[gene_id, 1]

                num_snps += 1

            df_ref_mt[:, gene_id] = np.divide(gene_data, num_snps)

            gene_mt = config.params_dict['genes_list'][0][genes_ids[gene_id]]
            names_mt.append(gene_mt)

        if int(config.params_dict['run_timer']) == 1:
            print('Time for data frame creating: ' + str(time.process_time() - start_df))

        df_ref_mt = df_ref_mt[:, : len(names_mt)]

        data_classes = []
        for item in target_samples_names_mt:
            if item in config.pop_person_dict[target_pop]:
                data_classes.append(target_pop)
            elif item in config.pop_person_dict[reference_pop]:
                data_classes.append(reference_pop)

        factor = pd.factorize(data_classes)
        y = factor[0]

        if int(config.params_dict['run_timer']) == 1:
            start_rf = time.process_time()

        clf = RandomForestClassifier(n_estimators=100)
        output = cross_validate(clf, df_ref_mt, y, cv=10, scoring='accuracy', return_estimator=True)
        accuracy = np.mean(output['test_score'])
        if int(config.params_dict['run_timer']) == 1:
            print('Total time for random forest ' + str(time.process_time() - start_rf))
            print('Mean accuracy: ' + str(accuracy))

        features_dict = dict((key, []) for key in names_mt)

        num_features = 0

        for idx, estimator in enumerate(output['estimator']):
            feature_importances = pd.DataFrame(estimator.feature_importances_,
                                               index=names_mt,
                                               columns=['importance']).sort_values('importance', ascending=False)

            features_names = list(feature_importances.index.values)
            features_values = list(feature_importances.values)
            for i in range(0, len(features_names)):
                features_dict[features_names[i]].append(features_values[i][0])

        for key in features_dict.keys():
            features_dict[key] = np.mean(features_dict[key])
            if features_dict[key] > 0.001:
                num_features += 1

        if accuracy >= float(config.params_dict['target_accuracy']):
            results.accuracy.append(accuracy)
            results.num_features.append(num_features)
            results.mt_genes.append(genes_ids)

        if int(config.params_dict['num_features']) > 0:
            features_dict = {k: v for k, v in sorted(features_dict.items(), reverse=True, key=lambda x: x[1])}
            results.features = features_dict


# type 3 - using mean frequencies over genes for sequential rf
def rf_type_3_mt(config, results):
    reference_part = float(config.params_dict['reference_part'])
    reference_pop = config.params_dict['reference_pop']
    target_pop = config.params_dict['target_pop']
    reference_size = int(len(config.pop_person_dict[reference_pop]) * reference_part)
    reference_list = random.sample(config.pop_person_dict[reference_pop], reference_size)

    genes_ids = config.params_dict['config_mt_genes'][0]
    genes_names = [config.params_dict['genes_list'][0][i] for i in genes_ids]

    print(';'.join(genes_names))

    if not hasattr(config, 'main_df'):

        # Reference group

        reference_frequencies = np.zeros((len(genes_ids), 2), dtype=np.float32)

        persons = config.person_index_dict

        target_samples_ids_mt = []
        target_samples_names_mt = []

        for sample_name in persons:
            if sample_name in reference_list:
                target_samples_ids_mt.append(persons[sample_name])
                target_samples_names_mt.append(sample_name)

        if int(config.params_dict['run_timer']) == 1:
            start_ref = time.process_time()

        for gene_id in range(0, len(genes_ids)):
            gene_index = config.data_position_dict[genes_names[gene_id]]
            for row in config.data[gene_index]:
                for i in range(0, len(target_samples_ids_mt)):
                    snp_mt = row[target_samples_ids_mt[i]]
                    if snp_mt == 0:
                        reference_frequencies[gene_id, 0] += 1
                    elif snp_mt == 1:
                        reference_frequencies[gene_id, 1] += 1

        for gene_id in range(0, len(genes_ids)):
            ref_sum = np.sum(reference_frequencies[gene_id, :])
            for i in range(0, 2):
                reference_frequencies[gene_id, i] = reference_frequencies[gene_id, i] / ref_sum

        if int(config.params_dict['run_timer']) == 1:
            print('Time for frequencies calculating: ' + str(time.process_time() - start_ref))

        # Remaining group

        target_samples_ids_mt = []
        target_samples_names_mt = []

        for sample_name in persons:
            if sample_name in config.pop_person_dict[target_pop] \
                    or sample_name in config.pop_person_dict[reference_pop]:
                target_samples_ids_mt.append(persons[sample_name])
                target_samples_names_mt.append(sample_name)

        num_frequencies = len(genes_ids)

        df_ref_mt = np.empty(shape=(len(target_samples_names_mt), num_frequencies), dtype=np.float32)

        if int(config.params_dict['run_timer']) == 1:
            start_df = time.process_time()

        names_mt = []
        gene_col_dict = {}

        for gene_id in range(0, len(genes_ids)):
            gene_index = config.data_position_dict[genes_names[gene_id]]
            gene_data = np.zeros(shape=len(target_samples_ids_mt), dtype=np.float32)
            num_snps = 0
            for row in config.data[gene_index]:
                if len(set(row[target_samples_ids_mt].tolist())) == 1:
                    continue
                for i in range(0, len(target_samples_ids_mt)):
                    snp_mt = row[target_samples_ids_mt[i]]
                    if snp_mt == 0:
                        gene_data[i] += 1 - reference_frequencies[gene_id, 0]
                    elif snp_mt == 1:
                        gene_data[i] += 1 - reference_frequencies[gene_id, 1]

                num_snps += 1

            df_ref_mt[:, gene_id] = np.divide(gene_data, num_snps)

            gene_mt = config.params_dict['genes_list'][0][genes_ids[gene_id]]
            names_mt.append(gene_mt)
            gene_col_dict[gene_mt] = gene_id

        if int(config.params_dict['run_timer']) == 1:
            print('Time for common data frame creating: ' + str(time.process_time() - start_df))

        data_classes = []
        for item in target_samples_names_mt:
            if item in config.pop_person_dict[target_pop]:
                data_classes.append(target_pop)
            elif item in config.pop_person_dict[reference_pop]:
                data_classes.append(reference_pop)

        factor = pd.factorize(data_classes)
        y = factor[0]

        config.main_df = df_ref_mt
        config.main_df_classes = y
        config.gene_col_dict = gene_col_dict
        save_df(config)

    else:
        df_ref_mt = config.main_df
        y = config.main_df_classes
        gene_col_dict = config.gene_col_dict
        names_mt = list(gene_col_dict.keys())

    if int(config.params_dict['run_timer']) == 1:
        start_rf = time.process_time()

    num_estimators = int(config.params_dict['num_estimators'])
    num_cv_runs = int(config.params_dict['num_cv_runs'])

    clf = RandomForestClassifier(n_estimators=num_estimators)
    output = cross_validate(clf, df_ref_mt, y, cv=num_cv_runs, scoring='accuracy', return_estimator=True)
    accuracy = np.mean(output['test_score'])
    if int(config.params_dict['run_timer']) == 1:
        print('Total time for random forest ' + str(time.process_time() - start_rf))

    features_dict = dict((key, []) for key in names_mt)

    num_features = 0

    for idx, estimator in enumerate(output['estimator']):
        feature_importances = pd.DataFrame(estimator.feature_importances_,
                                           index=names_mt,
                                           columns=['importance']).sort_values('importance', ascending=False)

        features_names = list(feature_importances.index.values)
        features_values = list(feature_importances.values)
        for i in range(0, len(features_names)):
            features_dict[features_names[i]].append(features_values[i][0])

    for key in features_dict.keys():
        features_dict[key] = np.mean(features_dict[key])
        if features_dict[key] > 0.001:
            num_features += 1

    if accuracy >= float(config.params_dict['target_accuracy']):
        results.accuracy.append(accuracy)
        results.num_features.append(num_features)
        results.mt_genes.append(genes_ids)

    features_dict = {k: v for k, v in sorted(features_dict.items(), reverse=True, key=lambda x: x[1])}
    results.features.append(features_dict)

    features_top = list(features_dict.keys())

    if config.params_dict['sequential_run_type'] == 'lin':
        features_counts = [i + 1 for i in range(0, int(config.params_dict['num_sequential_runs']))]
    elif config.params_dict['sequential_run_type'] == 'max':
        features_counts = [i + 1 for i in range(0, len(features_top) - 1)]
    else:
        features_counts = np.geomspace(1.0, len(features_top), int(config.params_dict['num_sequential_runs']),
                                       endpoint=True)
    features_counts = list(set([int(item) for item in features_counts]))
    features_counts.sort()

    for num_features in features_counts:
        if num_features % 10 == 0:
            print('Sequential random forest #' + str(num_features))

        curr_features = features_top[0:num_features]
        curr_features_ids = [gene_col_dict[feature] for feature in curr_features]
        curr_df = df_ref_mt[:, curr_features_ids].copy()

        clf = RandomForestClassifier(n_estimators=num_estimators)
        output = cross_validate(clf, curr_df, y, cv=num_cv_runs, scoring='accuracy', return_estimator=True)
        accuracy = np.mean(output['test_score'])

        if accuracy >= float(config.params_dict['target_accuracy']):
            results.accuracy.append(accuracy)
            results.mt_genes.append(curr_features_ids)

            features_dict = dict((key, []) for key in curr_features)

            for idx, estimator in enumerate(output['estimator']):
                feature_importances = pd.DataFrame(estimator.feature_importances_,
                                                   index=curr_features,
                                                   columns=['importance']).sort_values('importance', ascending=False)

                features_names = list(feature_importances.index.values)
                features_values = list(feature_importances.values)
                for i in range(0, len(features_names)):
                    features_dict[features_names[i]].append(features_values[i][0])

            for key in features_dict.keys():
                features_dict[key] = np.mean(features_dict[key])

            features_dict = {k: v for k, v in sorted(features_dict.items(), reverse=True, key=lambda x: x[1])}
            results.features.append(features_dict)


# type 4 - using genes variations for rf
def rf_type_4_mt(config, results):
    reference_pop = config.params_dict['reference_pop']
    target_pop = config.params_dict['target_pop']

    genes_ids_mt = config.params_dict['config_mt_genes'][0]
    genes_names_mt = [config.params_dict['genes_list'][0][i] for i in genes_ids_mt]

    print(';'.join(genes_names_mt))

    if not hasattr(config, 'main_df'):

        persons_mt = config.person_index_dict

        target_samples_ids_mt = []
        target_samples_names = []

        for sample_name in persons_mt:
            if sample_name in config.pop_person_dict[target_pop] \
                    or sample_name in config.pop_person_dict[reference_pop]:
                target_samples_ids_mt.append(persons_mt[sample_name])
                target_samples_names.append(sample_name)

        num_mt_snps = 0
        for gene_id_mt in range(0, len(genes_ids_mt)):
            gene_mt_index = config.data_position_dict[genes_names_mt[gene_id_mt]]
            for row_mt in config.data[gene_mt_index]:
                num_mt_snps += 1

        df_ref = np.empty(shape=(len(target_samples_names), num_mt_snps), dtype=np.float32)

        names = []
        snp_col_dict = {}

        if int(config.params_dict['run_timer']) == 1:
            start_df = time.process_time()

        snp_index = 0
        for gene_id_mt in range(0, len(genes_ids_mt)):
            gene_mt_index = config.data_position_dict[genes_names_mt[gene_id_mt]]
            print(config.params_dict['genes_list'][0][genes_ids_mt[gene_id_mt]])
            for row_mt_index in range(0, len(config.data[gene_mt_index])):
                row_mt = config.data[gene_mt_index][row_mt_index]
                if row_mt_index % 10 == 0:
                    print(str(row_mt_index))
                    if len(set(row_mt[target_samples_ids_mt].tolist())) == 1:
                        continue
                    for i in range(0, len(target_samples_ids_mt)):
                        snp_mt = row_mt[target_samples_ids_mt[i]]
                        if snp_mt == 0:
                            df_ref[i, snp_index] = 0
                        elif snp_mt == 1:
                            df_ref[i, snp_index] = 1

                    snp_index += 1

                    gene_mt = config.params_dict['genes_list'][0][genes_ids_mt[gene_id_mt]]
                    snp_mt = list(config.gene_snp_dict[gene_mt].keys())[row_mt_index]
                    name = snp_mt + ':' + gene_mt
                    names.append(name)
                    snp_col_dict[name] = snp_index

        df_ref = df_ref[:, :len(names)]

        if int(config.params_dict['run_timer']) == 1:
            print('Time for common data frame creating: ' + str(time.process_time() - start_df))

        data_classes = []
        for item in target_samples_names:
            if item in config.pop_person_dict[target_pop]:
                data_classes.append(target_pop)
            elif item in config.pop_person_dict[reference_pop]:
                data_classes.append(reference_pop)

        factor = pd.factorize(data_classes)
        y = factor[0]

        config.main_df = df_ref
        config.main_df_classes = y
        config.gene_col_dict = snp_col_dict
        save_df(config)

    else:
        df_ref = config.main_df
        y = config.main_df_classes
        gene_col_dict = config.gene_col_dict
        names = list(gene_col_dict.keys())

    if int(config.params_dict['run_timer']) == 1:
        start_rf = time.process_time()

    num_estimators = int(config.params_dict['num_estimators'])
    num_cv_runs = int(config.params_dict['num_cv_runs'])

    clf = RandomForestClassifier(n_estimators=num_estimators)
    output = cross_validate(clf, df_ref, y, cv=num_cv_runs, scoring='accuracy', return_estimator=True)
    accuracy = np.mean(output['test_score'])
    if int(config.params_dict['run_timer']) == 1:
        print('Total time for random forest ' + str(time.process_time() - start_rf))

    features_dict = dict((key, []) for key in names)

    num_features = 0

    for idx, estimator in enumerate(output['estimator']):
        feature_importances = pd.DataFrame(estimator.feature_importances_,
                                           index=names,
                                           columns=['importance']).sort_values('importance', ascending=False)

        features_names = list(feature_importances.index.values)
        features_values = list(feature_importances.values)
        for i in range(0, len(features_names)):
            features_dict[features_names[i]].append(features_values[i][0])

    for key in features_dict.keys():
        features_dict[key] = np.mean(features_dict[key])
        if features_dict[key] > 0:
            num_features += 1

    if accuracy >= float(config.params_dict['target_accuracy']):
        results.accuracy.append(accuracy)
        results.num_features.append(num_features)
        results.mt_genes.append(genes_ids_mt)
        results.nuc_genes.append([])

    features_dict = {k: v for k, v in sorted(features_dict.items(), reverse=True, key=lambda x: x[1])}
    results.features.append(features_dict)
