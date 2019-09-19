import itertools
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_validate

def unit_task(config, results):
    if config.params_dict['experiment_type'] == 'mt':
        task_mt(config, results)

def task_mt(config, results):
    reference_part = int(config.params_dict['reference_part'])
    reference_pop = config.params_dict['reference_pop']
    target_pop = config.params_dict['target_pop']
    reference_size = int(len(config.pop_person_dict[reference_pop]) * reference_part)
    reference_list = config.pop_person_dict[reference_pop][:reference_size]
    reference_frequencies = [0, 0]

    all_genes = config.params_dict['genes_ids_list']

    L = int(config.params_dict['k_mt'])

    for subset in itertools.combinations(all_genes, L):
        genes = list(subset)
        genes_ids = []
        for gene_id in range(0, len(genes)):
            genes_ids.append(config.data_position_dict[genes[gene_id]])

        print(';'.join(genes))

        # Reference group

        persons = config.person_index_dict

        target_samples_ids_mt = []
        target_samples_names_mt = []
        snp_samples_mt = {}

        for sample_name in persons:
            if sample_name in reference_list:
                target_samples_ids_mt.append(persons[sample_name])
                target_samples_names_mt.append(sample_name)
                snp_samples_mt[sample_name] = []

        number_mt_snps = 0

        for gene_id in genes_ids:
            for row in config.data[gene_id]:
                snp_data_mt = list(row[i] for i in target_samples_ids_mt)
                for id in range(0, len(snp_data_mt)):
                    if snp_data_mt[id] == '0':
                        reference_frequencies[0] += 1
                    elif snp_data_mt[id] == '1':
                        reference_frequencies[1] += 1

        reference_frequencies = [freq / sum(reference_frequencies) for freq in reference_frequencies]

        # Remaining group

        target_samples_ids_mt = []
        target_samples_names_mt = []

        for sample_name in persons:
            if sample_name in config.pop_person_dict[target_pop] \
                     or sample_name in config.pop_person_dict[reference_pop]:
                target_samples_ids_mt.append(persons[sample_name])
                target_samples_names_mt.append(sample_name)

        df_ref_mt = np.empty(shape=(len(target_samples_names_mt), number_mt_snps), dtype=float)
        names_mt = []

        line_count_mt = 0
        for gene_id in genes_ids:
            print('gene #' + str(gene_id) + ' processing')
            row_id = 0
            for row in config.data[gene_id]:
                snp_data_mt = list(row[i] for i in target_samples_ids_mt)

                combination_data = []

                for id in range(0, len(snp_data_mt)):
                    if snp_data_mt[id] == '0':
                        combination_data.append(1 - reference_frequencies[0])
                    elif snp_data_mt[id] == '1':
                        combination_data.append(1 - reference_frequencies[1])

                df_ref_mt[:, line_count_mt] = combination_data

                if len(set(combination_data)) > 1:
                    gene_mt = config.params_dict['genes_list'][gene_id]
                    snp_pos_mt = [index for name, index in config.gene_snp_dict[gene_mt].items() if index == row_id][0]
                    name_mt = gene_mt + '_' + snp_pos_mt
                    if name_mt not in names_mt:
                        names_mt.append(name_mt)

                line_count_mt += 1
                row_id += 1

        df_ref_mt = df_ref_mt[:, ~np.all(df_ref_mt[1:] == df_ref_mt[:-1], axis=0)]

        data_classes = []
        for item in target_samples_names_mt:
            for pop in target_pop:
                if item in config.pop_person_dict[target_pop] or item in config.pop_person_dict[reference_pop]:
                    data_classes.append(pop)

        factor = pd.factorize(data_classes)
        y = factor[0]

        clf = RandomForestClassifier(n_estimators=10)
        output = cross_validate(clf, df_ref_mt, y, cv=5, scoring='accuracy', return_estimator=True)
        accuracy = np.mean(output['test_score'])

        if accuracy >= int(config.params_dict['target_accuracy']):
            results.accuracy.append(accuracy)
            results.accuracy_mt_genes.append(genes)

        if int(config.params_dict['num_features']) > 0:

            features_dict = dict((key, []) for key in names_mt)

            for idx, estimator in enumerate(output['estimator']):
                feature_importances = pd.DataFrame(estimator.feature_importances_,
                                                   index=names_mt,
                                                   columns=['importance']).sort_values('importance', ascending=False)

                features_names = list(feature_importances.index.values)
                features_values = list(feature_importances.values)
                for id in range(0, len(features_names)):
                    features_dict[features_names[id]].append(features_values[0])

            for key in features_dict.keys():
                features_dict[key] = np.mean(features_dict[key])

            features_dict = {k: v for k, v in sorted(features_dict.items(), reverse=True, key=lambda x: x[1])}

            results.features = features_dict
