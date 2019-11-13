import pandas as pd
import numpy as np
import random
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_validate
import time

# type 0 - calculating of 6 reference frequencies
def rf_type_0_mt_nuc(config, results):
    reference_part = float(config.params_dict['reference_part'])
    reference_pop = config.params_dict['reference_pop']
    target_pop = config.params_dict['target_pop']
    reference_size = int(len(config.pop_person_dict[reference_pop]) * reference_part)
    reference_list = random.sample(config.pop_person_dict[reference_pop], reference_size)
    reference_frequencies = [0, 0, 0, 0, 0, 0]

    for task_id in range(0, len(config.params_dict['config_mt_genes'])):
        genes_ids_mt = config.params_dict['config_mt_genes'][task_id]
        genes_names_mt = [config.params_dict['genes_list'][0][i] for i in genes_ids_mt]

        genes_ids_nuc = config.params_dict['config_nuc_genes'][task_id]
        genes_names_nuc = [config.params_dict['genes_list'][1][i] for i in genes_ids_nuc]

        print(';'.join(genes_names_mt) + ';' + ';'.join(genes_names_nuc))

        # Reference group

        persons_mt = config.person_index_dict[0]
        persons_nuc = config.person_index_dict[1]

        target_samples_ids_mt = []
        target_samples_ids_nuc = []
        target_samples_names = []

        for sample_name in persons_mt:
            if sample_name in reference_list and sample_name in persons_nuc:
                target_samples_ids_mt.append(persons_mt[sample_name])
                target_samples_ids_nuc.append(persons_nuc[sample_name])
                target_samples_names.append(sample_name)

        number_snps_combinations = 0

        if int(config.params_dict['run_timer']) == 1:
            start_ref = time.process_time()

        for gene_id_mt in range(0, len(genes_ids_mt)):
            gene_mt_index = config.data_position_dict[genes_names_mt[gene_id_mt]]
            for gene_id_nuc in range(0, len(genes_ids_nuc)):
                gene_nuc_index = config.data_position_dict[genes_names_nuc[gene_id_nuc]]
                for row_mt in config.data[gene_mt_index]:
                    snp_data_mt = list(row_mt[i] for i in target_samples_ids_mt)
                    for row_nuc in config.data[gene_nuc_index]:
                        snp_data_nuc = list(row_nuc[i] for i in target_samples_ids_nuc)
                        for id in range(0, len(snp_data_mt)):
                            if snp_data_mt[id] == 0:
                                if snp_data_nuc[id] == 0:
                                    reference_frequencies[0] += 1
                                elif snp_data_nuc[id] == 1 or snp_data_nuc[id] == 2:
                                    reference_frequencies[1] += 1
                                elif snp_data_nuc[id] == 3:
                                    reference_frequencies[2] += 1
                            elif snp_data_mt[id] == 1:
                                if snp_data_nuc[id] == 0:
                                    reference_frequencies[3] += 1
                                elif snp_data_nuc[id] == 1 or snp_data_nuc[id] == 2:
                                    reference_frequencies[4] += 1
                                elif snp_data_nuc[id] == 3:
                                    reference_frequencies[5] += 1
                        number_snps_combinations += 1

        reference_frequencies = [freq / sum(reference_frequencies) for freq in reference_frequencies]

        if int(config.params_dict['run_timer']) == 1:
            print('Time for frequencies calculating: ' + str(time.process_time() - start_ref))

        # Remaining group

        target_samples_ids_mt = []
        target_samples_ids_nuc = []
        target_samples_names = []

        for sample_name in persons_mt:
            if sample_name in persons_nuc:
                if sample_name in config.pop_person_dict[target_pop] \
                        or sample_name in config.pop_person_dict[reference_pop]:
                    target_samples_ids_mt.append(persons_mt[sample_name])
                    target_samples_ids_nuc.append(persons_nuc[sample_name])
                    target_samples_names.append(sample_name)

        df_ref = np.empty(shape=(len(target_samples_names), number_snps_combinations), dtype=float)

        names = []
        line_count = 0

        if int(config.params_dict['run_timer']) == 1:
            start_df = time.process_time()

        for gene_id_mt in range(0, len(genes_ids_mt)):
            gene_mt_index = config.data_position_dict[genes_names_mt[gene_id_mt]]
            for gene_id_nuc in range(0, len(genes_ids_nuc)):
                gene_nuc_index = config.data_position_dict[genes_names_nuc[gene_id_nuc]]
                row_id_mt = 0
                for row_mt in config.data[gene_mt_index]:
                    snp_data_mt = list(row_mt[i] for i in target_samples_ids_mt)
                    row_id_nuc = 0
                    for row_nuc in config.data[gene_nuc_index]:
                        snp_data_nuc = list(row_nuc[i] for i in target_samples_ids_nuc)
                        if len(set(snp_data_mt)) == 1 and len(set(snp_data_nuc)) == 1:
                            continue
                        combination_data = []
                        for id in range(0, len(snp_data_mt)):
                            if snp_data_mt[id] == 0:
                                if snp_data_nuc[id] == 0:
                                    combination_data.append(1 - reference_frequencies[0])
                                elif snp_data_nuc[id] == 1 or snp_data_nuc[id] == 2:
                                    combination_data.append(1 - reference_frequencies[1])
                                elif snp_data_nuc[id] == 3:
                                    combination_data.append(1 - reference_frequencies[2])
                            if snp_data_mt[id] == 1:
                                if snp_data_nuc[id] == 0:
                                    combination_data.append(1 - reference_frequencies[3])
                                elif snp_data_nuc[id] == 1 or snp_data_nuc[id] == 2:
                                    combination_data.append(1 - reference_frequencies[4])
                                elif snp_data_nuc[id] == 3:
                                    combination_data.append(1 - reference_frequencies[5])

                        df_ref[:, line_count] = combination_data
                        if len(set(combination_data)) > 1:
                            gene_mt = config.params_dict['genes_list'][0][genes_ids_mt[gene_id_mt]]
                            snp_pos_mt = \
                                [name for name, index in config.gene_snp_dict[gene_mt].items() if index == row_id_mt][0]
                            gene_nuc = config.params_dict['genes_list'][1][genes_ids_nuc[gene_id_nuc]]
                            snp_pos_nuc = \
                                [name for name, index in config.gene_snp_dict[gene_nuc].items() if
                                 index == row_id_nuc][0]
                            name = gene_mt + '_' + snp_pos_mt + '_' + gene_nuc + '_' + snp_pos_nuc
                            if name not in names:
                                names.append(name)

                        line_count += 1
                        row_id_nuc += 1
                    row_id_mt += 1

        if int(config.params_dict['run_timer']) == 1:
            print('Time for data frame creating: ' + str(time.process_time() - start_df))
            print(df_ref.shape)

        df_ref = df_ref[:, : len(names)]

        if int(config.params_dict['run_timer']) == 1:
            print(df_ref.shape)

        data_classes = []
        for item in target_samples_names:
            if item in config.pop_person_dict[target_pop]:
                data_classes.append(target_pop)
            elif item in config.pop_person_dict[reference_pop]:
                data_classes.append(reference_pop)

        factor = pd.factorize(data_classes)
        y = factor[0]

        if int(config.params_dict['run_timer']) == 1:
            start_rf = time.process_time()

        clf = RandomForestClassifier(n_estimators=100)
        output = cross_validate(clf, df_ref, y, cv=10, scoring='accuracy', return_estimator=True)
        accuracy = np.mean(output['test_score'])

        if int(config.params_dict['run_timer']) == 1:
            score_time = output['score_time']
            print('Cross validation time: ' + ', '.join([str(item) for item in score_time]))
            print('Total time for random forest ' + str(time.process_time() - start_rf))
            print('Mean accuracy: ' + str(accuracy))

        features_dict = dict((key, []) for key in names)

        num_features = 0

        for idx, estimator in enumerate(output['estimator']):
            feature_importances = pd.DataFrame(estimator.feature_importances_,
                                               index=names,
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
            results.mt_genes.append(genes_ids_mt)
            results.nuc_genes.append(genes_ids_nuc)

        if int(config.params_dict['num_features']) > 0:
            features_dict = {k: v for k, v in sorted(features_dict.items(), reverse=True, key=lambda x: x[1])}
            results.features = features_dict

# type 1 - using genes variations for rf
def rf_type_1_mt_nuc(config, results):
    reference_pop = config.params_dict['reference_pop']
    target_pop = config.params_dict['target_pop']

    for task_id in range(0, len(config.params_dict['config_mt_genes'])):
        genes_ids_mt = config.params_dict['config_mt_genes'][task_id]
        genes_names_mt = [config.params_dict['genes_list'][0][i] for i in genes_ids_mt]

        genes_ids_nuc = config.params_dict['config_nuc_genes'][task_id]
        genes_names_nuc = [config.params_dict['genes_list'][1][i] for i in genes_ids_nuc]

        print(';'.join(genes_names_mt) + ';' + ';'.join(genes_names_nuc))

        persons_mt = config.person_index_dict[0]
        persons_nuc = config.person_index_dict[1]

        number_mt_snps = 0
        number_nuc_snps = 0
        for gene_id_mt in range(0, len(genes_ids_mt)):
            gene_mt_index = config.data_position_dict[genes_names_mt[gene_id_mt]]
            number_mt_snps += config.data[gene_mt_index].shape[0]
        for gene_id_nuc in range(0, len(genes_ids_nuc)):
            gene_nuc_index = config.data_position_dict[genes_names_nuc[gene_id_nuc]]
            number_nuc_snps += config.data[gene_nuc_index].shape[0]

        number_snps_combinations = number_mt_snps * number_nuc_snps

        target_samples_ids_mt = []
        target_samples_ids_nuc = []
        target_samples_names = []

        for sample_name in persons_mt:
            if sample_name in persons_nuc:
                if sample_name in config.pop_person_dict[target_pop] \
                        or sample_name in config.pop_person_dict[reference_pop]:
                    target_samples_ids_mt.append(persons_mt[sample_name])
                    target_samples_ids_nuc.append(persons_nuc[sample_name])
                    target_samples_names.append(sample_name)

        df_ref = np.empty(shape=(len(target_samples_names), number_snps_combinations), dtype=int)

        names = []
        line_count = 0

        if int(config.params_dict['run_timer']) == 1:
            start_df = time.process_time()

        for gene_id_mt in range(0, len(genes_ids_mt)):
            gene_mt_index = config.data_position_dict[genes_names_mt[gene_id_mt]]
            for gene_id_nuc in range(0, len(genes_ids_nuc)):
                gene_nuc_index = config.data_position_dict[genes_names_nuc[gene_id_nuc]]
                row_id_mt = 0
                for row_mt in config.data[gene_mt_index]:
                    snp_data_mt = list(row_mt[i] for i in target_samples_ids_mt)
                    row_id_nuc = 0
                    for row_nuc in config.data[gene_nuc_index]:
                        snp_data_nuc = list(row_nuc[i] for i in target_samples_ids_nuc)
                        if len(set(snp_data_mt)) == 1 and len(set(snp_data_nuc)) == 1:
                            continue
                        combination_data = []
                        for id in range(0, len(snp_data_mt)):
                            if snp_data_mt[id] == 0:
                                if snp_data_nuc[id] == 0:
                                    combination_data.append(0)
                                elif snp_data_nuc[id] == 1 or snp_data_nuc[id] == 2:
                                    combination_data.append(1)
                                elif snp_data_nuc[id] == 3:
                                    combination_data.append(2)
                            if snp_data_mt[id] == 1:
                                if snp_data_nuc[id] == 0:
                                    combination_data.append(3)
                                elif snp_data_nuc[id] == 1 or snp_data_nuc[id] == 2:
                                    combination_data.append(4)
                                elif snp_data_nuc[id] == 3:
                                    combination_data.append(5)

                        df_ref[:, line_count] = combination_data
                        if len(set(combination_data)) > 1:
                            gene_mt = config.params_dict['genes_list'][0][genes_ids_mt[gene_id_mt]]
                            snp_pos_mt = \
                                [name for name, index in config.gene_snp_dict[gene_mt].items() if index == row_id_mt][0]
                            gene_nuc = config.params_dict['genes_list'][1][genes_ids_nuc[gene_id_nuc]]
                            snp_pos_nuc = \
                                [name for name, index in config.gene_snp_dict[gene_nuc].items() if
                                 index == row_id_nuc][0]
                            name = gene_mt + '_' + snp_pos_mt + '_' + gene_nuc + '_' + snp_pos_nuc
                            if name not in names:
                                names.append(name)

                        line_count += 1
                        row_id_nuc += 1
                    row_id_mt += 1

        if int(config.params_dict['run_timer']) == 1:
            print('Time for data frame creating: ' + str(time.process_time() - start_df))
            print(df_ref.shape)

        df_ref = df_ref[:, : len(names)]

        if int(config.params_dict['run_timer']) == 1:
            print(df_ref.shape)

        data_classes = []
        for item in target_samples_names:
            if item in config.pop_person_dict[target_pop]:
                data_classes.append(target_pop)
            elif item in config.pop_person_dict[reference_pop]:
                data_classes.append(reference_pop)

        factor = pd.factorize(data_classes)
        y = factor[0]

        if int(config.params_dict['run_timer']) == 1:
            start_rf = time.process_time()

        clf = RandomForestClassifier(n_estimators=100)
        output = cross_validate(clf, df_ref, y, cv=10, scoring='accuracy', return_estimator=True)
        accuracy = np.mean(output['test_score'])

        if int(config.params_dict['run_timer']) == 1:
            score_time = output['score_time']
            print('Cross validation time: ' + ', '.join([str(item) for item in score_time]))
            print('Total time for random forest: ' + str(time.process_time() - start_rf))
            print('Mean accuracy: ' + str(accuracy))

        features_dict = dict((key, []) for key in names)

        num_features = 0

        for idx, estimator in enumerate(output['estimator']):
            feature_importances = pd.DataFrame(estimator.feature_importances_,
                                               index=names,
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
            results.mt_genes.append(genes_ids_mt)
            results.nuc_genes.append(genes_ids_nuc)

        if int(config.params_dict['num_features']) > 0:
            features_dict = {k: v for k, v in sorted(features_dict.items(), reverse=True, key=lambda x: x[1])}
            results.features = features_dict


# type 2 - using mean frequencies over genes for rf
def rf_type_2_mt_nuc(config, results):
    reference_part = float(config.params_dict['reference_part'])
    reference_pop = config.params_dict['reference_pop']
    target_pop = config.params_dict['target_pop']
    reference_size = int(len(config.pop_person_dict[reference_pop]) * reference_part)
    reference_list = random.sample(config.pop_person_dict[reference_pop], reference_size)

    for task_id in range(0, len(config.params_dict['config_mt_genes'])):
        genes_ids_mt = config.params_dict['config_mt_genes'][task_id]
        genes_names_mt = [config.params_dict['genes_list'][0][i] for i in genes_ids_mt]

        genes_ids_nuc = config.params_dict['config_nuc_genes'][task_id]
        genes_names_nuc = [config.params_dict['genes_list'][1][i] for i in genes_ids_nuc]

        print(';'.join(genes_names_mt) + ';' + ';'.join(genes_names_nuc))

        # Reference group

        reference_frequencies = np.zeros((len(genes_ids_mt) * len(genes_ids_nuc), 6), dtype=np.float32)

        persons_mt = config.person_index_dict[0]
        persons_nuc = config.person_index_dict[1]

        target_samples_ids_mt = []
        target_samples_ids_nuc = []
        target_samples_names = []

        for sample_name in persons_mt:
            if sample_name in reference_list and sample_name in persons_nuc:
                target_samples_ids_mt.append(persons_mt[sample_name])
                target_samples_ids_nuc.append(persons_nuc[sample_name])
                target_samples_names.append(sample_name)

        if int(config.params_dict['run_timer']) == 1:
            start_ref = time.process_time()

        for gene_id_mt in range(0, len(genes_ids_mt)):
            gene_mt_index = config.data_position_dict[genes_names_mt[gene_id_mt]]
            for gene_id_nuc in range(0, len(genes_ids_nuc)):
                gene_nuc_index = config.data_position_dict[genes_names_nuc[gene_id_nuc]]
                common_index = gene_id_mt * len(genes_ids_nuc) + gene_id_nuc
                for row_mt in config.data[gene_mt_index]:
                    snp_data_mt = list(row_mt[i] for i in target_samples_ids_mt)
                    for row_nuc in config.data[gene_nuc_index]:
                        snp_data_nuc = list(row_nuc[i] for i in target_samples_ids_nuc)
                        for id in range(0, len(snp_data_mt)):
                            if snp_data_mt[id] == 0:
                                if snp_data_nuc[id] == 0:
                                    reference_frequencies[common_index, 0] += 1
                                elif snp_data_nuc[id] == 1 or snp_data_nuc[id] == 2:
                                    reference_frequencies[common_index, 1] += 1
                                elif snp_data_nuc[id] == 3:
                                    reference_frequencies[common_index, 2] += 1
                            elif snp_data_mt[id] == 1:
                                if snp_data_nuc[id] == 0:
                                    reference_frequencies[common_index, 3] += 1
                                elif snp_data_nuc[id] == 1 or snp_data_nuc[id] == 2:
                                    reference_frequencies[common_index, 4] += 1
                                elif snp_data_nuc[id] == 3:
                                    reference_frequencies[common_index, 5] += 1

        for gene_id in range(0, len(genes_ids_mt) * len(genes_ids_nuc)):
            ref_sum = np.sum(reference_frequencies[gene_id, :])
            for i in range(0, 6):
                reference_frequencies[gene_id, i] = reference_frequencies[gene_id, i] / ref_sum

        if int(config.params_dict['run_timer']) == 1:
            print('Time for frequencies calculating: ' + str(time.process_time() - start_ref))

        # Remaining group

        target_samples_ids_mt = []
        target_samples_ids_nuc = []
        target_samples_names = []

        for sample_name in persons_mt:
            if sample_name in persons_nuc:
                if sample_name in config.pop_person_dict[target_pop] \
                        or sample_name in config.pop_person_dict[reference_pop]:
                    target_samples_ids_mt.append(persons_mt[sample_name])
                    target_samples_ids_nuc.append(persons_nuc[sample_name])
                    target_samples_names.append(sample_name)

        df_ref = np.empty(shape=(len(target_samples_names), number_snps_combinations), dtype=float)

        names = []
        line_count = 0

        if int(config.params_dict['run_timer']) == 1:
            start_df = time.process_time()

        for gene_id_mt in range(0, len(genes_ids_mt)):
            gene_mt_index = config.data_position_dict[genes_names_mt[gene_id_mt]]
            for gene_id_nuc in range(0, len(genes_ids_nuc)):
                gene_nuc_index = config.data_position_dict[genes_names_nuc[gene_id_nuc]]
                row_id_mt = 0
                for row_mt in config.data[gene_mt_index]:
                    snp_data_mt = list(row_mt[i] for i in target_samples_ids_mt)
                    row_id_nuc = 0
                    for row_nuc in config.data[gene_nuc_index]:
                        snp_data_nuc = list(row_nuc[i] for i in target_samples_ids_nuc)
                        if len(set(snp_data_mt)) == 1 and len(set(snp_data_nuc)) == 1:
                            continue
                        combination_data = []
                        for id in range(0, len(snp_data_mt)):
                            if snp_data_mt[id] == 0:
                                if snp_data_nuc[id] == 0:
                                    combination_data.append(1 - reference_frequencies[0])
                                elif snp_data_nuc[id] == 1 or snp_data_nuc[id] == 2:
                                    combination_data.append(1 - reference_frequencies[1])
                                elif snp_data_nuc[id] == 3:
                                    combination_data.append(1 - reference_frequencies[2])
                            if snp_data_mt[id] == 1:
                                if snp_data_nuc[id] == 0:
                                    combination_data.append(1 - reference_frequencies[3])
                                elif snp_data_nuc[id] == 1 or snp_data_nuc[id] == 2:
                                    combination_data.append(1 - reference_frequencies[4])
                                elif snp_data_nuc[id] == 3:
                                    combination_data.append(1 - reference_frequencies[5])

                        df_ref[:, line_count] = combination_data
                        if len(set(combination_data)) > 1:
                            gene_mt = config.params_dict['genes_list'][0][genes_ids_mt[gene_id_mt]]
                            snp_pos_mt = \
                                [name for name, index in config.gene_snp_dict[gene_mt].items() if index == row_id_mt][0]
                            gene_nuc = config.params_dict['genes_list'][1][genes_ids_nuc[gene_id_nuc]]
                            snp_pos_nuc = \
                                [name for name, index in config.gene_snp_dict[gene_nuc].items() if
                                 index == row_id_nuc][0]
                            name = gene_mt + '_' + snp_pos_mt + '_' + gene_nuc + '_' + snp_pos_nuc
                            if name not in names:
                                names.append(name)

                        line_count += 1
                        row_id_nuc += 1
                    row_id_mt += 1

        if int(config.params_dict['run_timer']) == 1:
            print('Time for data frame creating: ' + str(time.process_time() - start_df))
            print(df_ref.shape)

        df_ref = df_ref[:, : len(names)]

        if int(config.params_dict['run_timer']) == 1:
            print(df_ref.shape)

        data_classes = []
        for item in target_samples_names:
            if item in config.pop_person_dict[target_pop]:
                data_classes.append(target_pop)
            elif item in config.pop_person_dict[reference_pop]:
                data_classes.append(reference_pop)

        factor = pd.factorize(data_classes)
        y = factor[0]

        if int(config.params_dict['run_timer']) == 1:
            start_rf = time.process_time()

        clf = RandomForestClassifier(n_estimators=100)
        output = cross_validate(clf, df_ref, y, cv=10, scoring='accuracy', return_estimator=True)
        accuracy = np.mean(output['test_score'])

        if int(config.params_dict['run_timer']) == 1:
            score_time = output['score_time']
            print('Cross validation time: ' + ', '.join([str(item) for item in score_time]))
            print('Total time for random forest ' + str(time.process_time() - start_rf))
            print('Mean accuracy: ' + str(accuracy))

        features_dict = dict((key, []) for key in names)

        num_features = 0

        for idx, estimator in enumerate(output['estimator']):
            feature_importances = pd.DataFrame(estimator.feature_importances_,
                                               index=names,
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
            results.mt_genes.append(genes_ids_mt)
            results.nuc_genes.append(genes_ids_nuc)

        if int(config.params_dict['num_features']) > 0:
            features_dict = {k: v for k, v in sorted(features_dict.items(), reverse=True, key=lambda x: x[1])}
            results.features = features_dict
