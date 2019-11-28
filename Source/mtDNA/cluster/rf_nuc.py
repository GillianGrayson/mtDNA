import pandas as pd
import numpy as np
import random
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_validate
import time

# type 0 - calculating of 3 reference frequencies
def rf_type_0_nuc(config, results):
    reference_part = float(config.params_dict['reference_part'])
    reference_pop = config.params_dict['reference_pop']
    target_pop = config.params_dict['target_pop']
    reference_size = int(len(config.pop_person_dict[reference_pop]) * reference_part)
    reference_list = random.sample(config.pop_person_dict[reference_pop], reference_size)
    reference_frequencies = [0, 0, 0]

    for task in config.params_dict['config_nuc_genes']:

        genes_ids = list(task)
        genes_names = [config.params_dict['genes_list'][0][i] for i in genes_ids]

        print(';'.join(genes_names))

        # Reference group

        persons = config.person_index_dict

        target_samples_ids_nuc = []
        target_samples_names_nuc = []

        for sample_name in persons:
            if sample_name in reference_list:
                target_samples_ids_nuc.append(persons[sample_name])
                target_samples_names_nuc.append(sample_name)

        number_nuc_snps = 0

        if int(config.params_dict['run_timer']) == 1:
            start_ref = time.process_time()

        for gene_id in range(0, len(genes_ids)):
            gene_index = config.data_position_dict[genes_names[gene_id]]
            for row in config.data[gene_index]:
                snp_data_nuc = list(row[i] for i in target_samples_ids_nuc)
                for id in range(0, len(snp_data_nuc)):
                    if snp_data_nuc[id] == 0:
                        reference_frequencies[0] += 1
                    elif snp_data_nuc[id] == 1 or snp_data_nuc[id] == 2:
                        reference_frequencies[1] += 1
                    elif snp_data_nuc[id] == 3:
                        reference_frequencies[2] += 1
                number_nuc_snps += 1

        reference_frequencies = [freq / sum(reference_frequencies) for freq in reference_frequencies]

        if int(config.params_dict['run_timer']) == 1:
            print('Time for frequencies calculating: ' + str(time.process_time() - start_ref))

        # Remaining group

        target_samples_ids_nuc = []
        target_samples_names_nuc = []

        for sample_name in persons:
            if sample_name in config.pop_person_dict[target_pop] \
                    or sample_name in config.pop_person_dict[reference_pop]:
                target_samples_ids_nuc.append(persons[sample_name])
                target_samples_names_nuc.append(sample_name)

        df_ref_nuc = np.empty(shape=(len(target_samples_names_nuc), number_nuc_snps), dtype=float)
        names_nuc = []

        line_count_nuc = 0

        if int(config.params_dict['run_timer']) == 1:
            start_df = time.process_time()

        for gene_id in range(0, len(genes_ids)):
            gene_index = config.data_position_dict[genes_names[gene_id]]
            row_id = 0
            for row in config.data[gene_index]:
                snp_data_nuc = list(row[i] for i in target_samples_ids_nuc)
                if len(set(snp_data_nuc)) == 1:
                    continue

                combination_data = []

                for id in range(0, len(snp_data_nuc)):
                    if snp_data_nuc[id] == 0:
                        combination_data.append(1 - reference_frequencies[0])
                    elif snp_data_nuc[id] == 1 or snp_data_nuc[id] == 2:
                        combination_data.append(1 - reference_frequencies[1])
                    elif snp_data_nuc[id] == 3:
                        combination_data.append(1 - reference_frequencies[2])

                df_ref_nuc[:, line_count_nuc] = combination_data

                if len(set(combination_data)) > 1:
                    gene_nuc = config.params_dict['genes_list'][0][genes_ids[gene_id]]
                    snp_pos_nuc = [name for name, index in config.gene_snp_dict[gene_nuc].items() if index == row_id][0]
                    name_nuc = gene_nuc + '_' + snp_pos_nuc
                    if name_nuc not in names_nuc:
                        names_nuc.append(name_nuc)

                line_count_nuc += 1
                row_id += 1

        if int(config.params_dict['run_timer']) == 1:
            print('Time for data frame creating: ' + str(time.process_time() - start_df))

        df_ref_nuc = df_ref_nuc[:, : len(names_nuc)]

        data_classes = []
        for item in target_samples_names_nuc:
            if item in config.pop_person_dict[target_pop]:
                data_classes.append(target_pop)
            elif item in config.pop_person_dict[reference_pop]:
                data_classes.append(reference_pop)

        factor = pd.factorize(data_classes)
        y = factor[0]

        if int(config.params_dict['run_timer']) == 1:
            start_rf = time.process_time()

        clf = RandomForestClassifier(n_estimators=100)
        output = cross_validate(clf, df_ref_nuc, y, cv=10, scoring='accuracy', return_estimator=True)
        accuracy = np.mean(output['test_score'])

        if int(config.params_dict['run_timer']) == 1:
            score_time = output['score_time']
            print('Cross validation time: ' + ', '.join([str(item) for item in score_time]))
            print('Total time for random forest ' + str(time.process_time() - start_rf))
            print('Mean accuracy: ' + str(accuracy))

        features_dict = dict((key, []) for key in names_nuc)

        num_features = 0

        for idx, estimator in enumerate(output['estimator']):
            feature_importances = pd.DataFrame(estimator.feature_importances_,
                                               index=names_nuc,
                                               columns=['importance']).sort_values('importance', ascending=False)

            features_names = list(feature_importances.index.values)
            features_values = list(feature_importances.values)
            for id in range(0, len(features_names)):
                features_dict[features_names[id]].append(features_values[id][0])

        for key in features_dict.keys():
            features_dict[key] = np.mean(features_dict[key])
            if features_dict[key] > 0:
                num_features += 1

        if accuracy >= float(config.params_dict['target_accuracy']):
            results.accuracy.append(accuracy)
            results.num_features.append(num_features)
            results.nuc_genes.append(genes_ids)

        if int(config.params_dict['num_features']) > 0:
            features_dict = {k: v for k, v in sorted(features_dict.items(), reverse=True, key=lambda x: x[1])}
            results.features = features_dict

# type 1 - using genes variations for rf
def rf_type_1_nuc(config, results):
    reference_pop = config.params_dict['reference_pop']
    target_pop = config.params_dict['target_pop']

    for task in config.params_dict['config_nuc_genes']:

        genes_ids = list(task)
        genes_names = [config.params_dict['genes_list'][0][i] for i in genes_ids]

        print(';'.join(genes_names))

        persons = config.person_index_dict

        number_nuc_snps = 0
        for gene_id in range(0, len(genes_ids)):
            gene_index = config.data_position_dict[genes_names[gene_id]]
            number_nuc_snps += config.data[gene_index].shape[0]

        target_samples_ids_nuc = []
        target_samples_names_nuc = []

        for sample_name in persons:
            if sample_name in config.pop_person_dict[target_pop] \
                    or sample_name in config.pop_person_dict[reference_pop]:
                target_samples_ids_nuc.append(persons[sample_name])
                target_samples_names_nuc.append(sample_name)

        df_ref_nuc = np.empty(shape=(len(target_samples_names_nuc), number_nuc_snps), dtype=int)
        names_nuc = []

        line_count_nuc = 0

        if int(config.params_dict['run_timer']) == 1:
            start_df = time.process_time()

        for gene_id in range(0, len(genes_ids)):
            gene_index = config.data_position_dict[genes_names[gene_id]]
            row_id = 0
            for row in config.data[gene_index]:
                snp_data_nuc = list(row[i] for i in target_samples_ids_nuc)
                if len(set(snp_data_nuc)) == 1:
                    continue

                combination_data = []

                for id in range(0, len(snp_data_nuc)):
                    if snp_data_nuc[id] == 0:
                        combination_data.append(0)
                    elif snp_data_nuc[id] == 1 or snp_data_nuc[id] == 2:
                        combination_data.append(1)
                    elif snp_data_nuc[id] == 3:
                        combination_data.append(2)

                df_ref_nuc[:, line_count_nuc] = combination_data

                if len(set(combination_data)) > 1:
                    gene_nuc = config.params_dict['genes_list'][0][genes_ids[gene_id]]
                    snp_pos_nuc = [name for name, index in config.gene_snp_dict[gene_nuc].items() if index == row_id][0]
                    name_nuc = gene_nuc + '_' + snp_pos_nuc
                    if name_nuc not in names_nuc:
                        names_nuc.append(name_nuc)

                line_count_nuc += 1
                row_id += 1

        if int(config.params_dict['run_timer']) == 1:
            print('Time for data frame creating: ' + str(time.process_time() - start_df))

        df_ref_nuc = df_ref_nuc[:, : len(names_nuc)]

        data_classes = []
        for item in target_samples_names_nuc:
            if item in config.pop_person_dict[target_pop]:
                data_classes.append(target_pop)
            elif item in config.pop_person_dict[reference_pop]:
                data_classes.append(reference_pop)

        factor = pd.factorize(data_classes)
        y = factor[0]

        if int(config.params_dict['run_timer']) == 1:
            start_rf = time.process_time()

        clf = RandomForestClassifier(n_estimators=100)
        output = cross_validate(clf, df_ref_nuc, y, cv=10, scoring='accuracy', return_estimator=True)
        accuracy = np.mean(output['test_score'])

        if int(config.params_dict['run_timer']) == 1:
            score_time = output['score_time']
            print('Cross validation time: ' + ', '.join([str(item) for item in score_time]))
            print('Total time for random forest ' + str(time.process_time() - start_rf))
            print('Mean accuracy: ' + str(accuracy))

        features_dict = dict((key, []) for key in names_nuc)

        num_features = 0

        for idx, estimator in enumerate(output['estimator']):
            feature_importances = pd.DataFrame(estimator.feature_importances_,
                                               index=names_nuc,
                                               columns=['importance']).sort_values('importance', ascending=False)

            features_names = list(feature_importances.index.values)
            features_values = list(feature_importances.values)
            for id in range(0, len(features_names)):
                features_dict[features_names[id]].append(features_values[id][0])

        for key in features_dict.keys():
            features_dict[key] = np.mean(features_dict[key])
            if features_dict[key] > 0:
                num_features += 1

        if accuracy >= float(config.params_dict['target_accuracy']):
            results.accuracy.append(accuracy)
            results.num_features.append(num_features)
            results.nuc_genes.append(genes_ids)

        if int(config.params_dict['num_features']) > 0:
            features_dict = {k: v for k, v in sorted(features_dict.items(), reverse=True, key=lambda x: x[1])}
            results.features = features_dict


# type 2 - using mean frequencies over genes for rf
def rf_type_2_nuc(config, results):
    reference_part = float(config.params_dict['reference_part'])
    reference_pop = config.params_dict['reference_pop']
    target_pop = config.params_dict['target_pop']
    reference_size = int(len(config.pop_person_dict[reference_pop]) * reference_part)
    reference_list = random.sample(config.pop_person_dict[reference_pop], reference_size)

    for task in config.params_dict['config_nuc_genes']:

        genes_ids = list(task)
        genes_names = [config.params_dict['genes_list'][0][i] for i in genes_ids]

        print(';'.join(genes_names))

        # Reference group

        reference_frequencies = np.zeros((len(genes_ids), 3), dtype=np.float32)

        persons = config.person_index_dict

        target_samples_ids_nuc = []
        target_samples_names_nuc = []

        for sample_name in persons:
            if sample_name in reference_list:
                target_samples_ids_nuc.append(persons[sample_name])
                target_samples_names_nuc.append(sample_name)

        if int(config.params_dict['run_timer']) == 1:
            start_ref = time.process_time()

        for gene_id in range(0, len(genes_ids)):
            gene_index = config.data_position_dict[genes_names[gene_id]]
            for row in config.data[gene_index]:
                snp_data_nuc = list(row[i] for i in target_samples_ids_nuc)
                for i in range(0, len(snp_data_nuc)):
                    if snp_data_nuc[i] == 0:
                        reference_frequencies[gene_id, 0] += 1
                    elif snp_data_nuc[i] == 1 or snp_data_nuc[i] == 2:
                        reference_frequencies[gene_id, 1] += 1
                    elif snp_data_nuc[i] == 3:
                        reference_frequencies[gene_id, 2] += 1

        for gene_id in range(0, len(genes_ids)):
            ref_sum = np.sum(reference_frequencies[gene_id, :])
            for i in range(0, 3):
                reference_frequencies[gene_id, i] = reference_frequencies[gene_id, i] / ref_sum

        if int(config.params_dict['run_timer']) == 1:
            print('Time for frequencies calculating: ' + str(time.process_time() - start_ref))

        # Remaining group

        target_samples_ids_nuc = []
        target_samples_names_nuc = []

        for sample_name in persons:
            if sample_name in config.pop_person_dict[target_pop] \
                    or sample_name in config.pop_person_dict[reference_pop]:
                target_samples_ids_nuc.append(persons[sample_name])
                target_samples_names_nuc.append(sample_name)

        num_frequencies = len(genes_ids)

        df_ref_nuc = np.empty(shape=(len(target_samples_names_nuc), num_frequencies), dtype=np.float32)

        if int(config.params_dict['run_timer']) == 1:
            start_df = time.process_time()

        names_nuc = []

        for gene_id in range(0, len(genes_ids)):
            gene_index = config.data_position_dict[genes_names[gene_id]]
            gene_data = np.zeros(shape=len(target_samples_ids_nuc), dtype=np.float32)
            num_snps = 0
            for row in config.data[gene_index]:
                snp_data_nuc = list(row[i] for i in target_samples_ids_nuc)
                if len(set(snp_data_nuc)) == 1:
                    continue

                for i in range(0, len(snp_data_nuc)):
                    if snp_data_nuc[i] == 0:
                        gene_data[i] += 1 - reference_frequencies[gene_id, 0]
                    elif snp_data_nuc[i] == 1 or snp_data_nuc[i] == 2:
                        gene_data[i] += 1 - reference_frequencies[gene_id, 1]
                    elif snp_data_nuc[i] == 3:
                        gene_data[i] += 1 - reference_frequencies[gene_id, 2]

                num_snps += 1

            df_ref_nuc[:, gene_id] = np.divide(gene_data, num_snps)

            gene_nuc = config.params_dict['genes_list'][0][genes_ids[gene_id]]
            names_nuc.append(gene_nuc)

        if int(config.params_dict['run_timer']) == 1:
            print('Time for data frame creating: ' + str(time.process_time() - start_df))

        df_ref_nuc = df_ref_nuc[:, : len(names_nuc)]

        data_classes = []
        for item in target_samples_names_nuc:
            if item in config.pop_person_dict[target_pop]:
                data_classes.append(target_pop)
            elif item in config.pop_person_dict[reference_pop]:
                data_classes.append(reference_pop)

        factor = pd.factorize(data_classes)
        y = factor[0]

        if int(config.params_dict['run_timer']) == 1:
            start_rf = time.process_time()

        clf = RandomForestClassifier(n_estimators=100)
        output = cross_validate(clf, df_ref_nuc, y, cv=10, scoring='accuracy', return_estimator=True)
        accuracy = np.mean(output['test_score'])
        if int(config.params_dict['run_timer']) == 1:
            print('Total time for random forest ' + str(time.process_time() - start_rf))
            print('Mean accuracy: ' + str(accuracy))

        features_dict = dict((key, []) for key in names_nuc)

        num_features = 0

        for idx, estimator in enumerate(output['estimator']):
            feature_importances = pd.DataFrame(estimator.feature_importances_,
                                               index=names_nuc,
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
            results.nuc_genes.append(genes_ids)

        if int(config.params_dict['num_features']) > 0:
            features_dict = {k: v for k, v in sorted(features_dict.items(), reverse=True, key=lambda x: x[1])}
            results.features = features_dict


# type 3 - using mean frequencies over genes for sequential rf
def rf_type_3_nuc(config, results):
    reference_part = float(config.params_dict['reference_part'])
    reference_pop = config.params_dict['reference_pop']
    target_pop = config.params_dict['target_pop']
    reference_size = int(len(config.pop_person_dict[reference_pop]) * reference_part)
    reference_list = random.sample(config.pop_person_dict[reference_pop], reference_size)

    genes_ids = config.params_dict['config_nuc_genes'][0]
    genes_names = [config.params_dict['genes_list'][0][i] for i in genes_ids]

    print(';'.join(genes_names))

    # Reference group

    reference_frequencies = np.zeros((len(genes_ids), 3), dtype=np.float32)

    persons = config.person_index_dict

    target_samples_ids_nuc = []
    target_samples_names_nuc = []

    for sample_name in persons:
        if sample_name in reference_list:
            target_samples_ids_nuc.append(persons[sample_name])
            target_samples_names_nuc.append(sample_name)

    if int(config.params_dict['run_timer']) == 1:
        start_ref = time.process_time()

    for gene_id in range(0, len(genes_ids)):
        gene_index = config.data_position_dict[genes_names[gene_id]]
        for row in config.data[gene_index]:
            snp_data_nuc = list(row[i] for i in target_samples_ids_nuc)
            for i in range(0, len(snp_data_nuc)):
                if snp_data_nuc[i] == 0:
                    reference_frequencies[gene_id, 0] += 1
                elif snp_data_nuc[i] == 1 or snp_data_nuc[i] == 2:
                    reference_frequencies[gene_id, 1] += 1
                elif snp_data_nuc[i] == 3:
                    reference_frequencies[gene_id, 2] += 1

    for gene_id in range(0, len(genes_ids)):
        ref_sum = np.sum(reference_frequencies[gene_id, :])
        for i in range(0, 3):
            reference_frequencies[gene_id, i] = reference_frequencies[gene_id, i] / ref_sum

    if int(config.params_dict['run_timer']) == 1:
        print('Time for frequencies calculating: ' + str(time.process_time() - start_ref))

    # Remaining group

    target_samples_ids_nuc = []
    target_samples_names_nuc = []

    for sample_name in persons:
        if sample_name in config.pop_person_dict[target_pop] \
                or sample_name in config.pop_person_dict[reference_pop]:
            target_samples_ids_nuc.append(persons[sample_name])
            target_samples_names_nuc.append(sample_name)

    num_frequencies = len(genes_ids)

    df_ref_nuc = np.empty(shape=(len(target_samples_names_nuc), num_frequencies), dtype=np.float32)

    if int(config.params_dict['run_timer']) == 1:
        start_df = time.process_time()

    names_nuc = []
    gene_col_dict = {}

    for gene_id in range(0, len(genes_ids)):
        gene_index = config.data_position_dict[genes_names[gene_id]]
        gene_data = np.zeros(shape=len(target_samples_ids_nuc), dtype=np.float32)
        num_snps = 0
        for row in config.data[gene_index]:
            snp_data_nuc = list(row[i] for i in target_samples_ids_nuc)
            if len(set(snp_data_nuc)) == 1:
                continue

            for i in range(0, len(snp_data_nuc)):
                if snp_data_nuc[i] == 0:
                    gene_data[i] += 1 - reference_frequencies[gene_id, 0]
                elif snp_data_nuc[i] == 1 or snp_data_nuc[i] == 2:
                    gene_data[i] += 1 - reference_frequencies[gene_id, 1]
                elif snp_data_nuc[i] == 3:
                    gene_data[i] += 1 - reference_frequencies[gene_id, 2]

            num_snps += 1

        df_ref_nuc[:, gene_id] = np.divide(gene_data, num_snps)

        gene_nuc = config.params_dict['genes_list'][0][genes_ids[gene_id]]
        names_nuc.append(gene_nuc)
        gene_col_dict[gene_nuc] = gene_id

    if int(config.params_dict['run_timer']) == 1:
        print('Time for common data frame creating: ' + str(time.process_time() - start_df))

    df_ref_nuc = df_ref_nuc[:, : len(names_nuc)]

    data_classes = []
    for item in target_samples_names_nuc:
        if item in config.pop_person_dict[target_pop]:
            data_classes.append(target_pop)
        elif item in config.pop_person_dict[reference_pop]:
            data_classes.append(reference_pop)

    factor = pd.factorize(data_classes)
    y = factor[0]

    if int(config.params_dict['run_timer']) == 1:
        start_rf = time.process_time()

    config.main_df = df_ref_nuc
    config.main_df_classes = y

    clf = RandomForestClassifier(n_estimators=100)
    output = cross_validate(clf, df_ref_nuc, y, cv=10, scoring='accuracy', return_estimator=True)
    accuracy = np.mean(output['test_score'])
    if int(config.params_dict['run_timer']) == 1:
        print('Total time for random forest ' + str(time.process_time() - start_rf))

    features_dict = dict((key, []) for key in names_nuc)

    num_features = 0

    for idx, estimator in enumerate(output['estimator']):
        feature_importances = pd.DataFrame(estimator.feature_importances_,
                                           index=names_nuc,
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
        results.nuc_genes.append(genes_ids)

    features_dict = {k: v for k, v in sorted(features_dict.items(), reverse=True, key=lambda x: x[1])}
    results.features = features_dict

    features_top = list(features_dict.keys())

    for num_features in range(1, len(features_top)):
        if num_features % 10 == 0:
            print('Sequential random forest #' + str(num_features))

        curr_features = features_top[0:num_features]
        curr_features_ids = [gene_col_dict[feature] for feature in curr_features]
        curr_df = df_ref_nuc[:, curr_features_ids].copy()

        clf = RandomForestClassifier(n_estimators=100)
        output = cross_validate(clf, curr_df, y, cv=10, scoring='accuracy', return_estimator=True)
        accuracy = np.mean(output['test_score'])

        if accuracy >= float(config.params_dict['target_accuracy']):
            results.accuracy.append(accuracy)
            results.nuc_genes.append(curr_features_ids)
