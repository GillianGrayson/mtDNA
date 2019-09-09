import os
import itertools
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_validate


def rf_mt(config_dict):
    data_path = config_dict['data_path'][0]
    result_path = config_dict['result_path'][0]
    ref_pop = config_dict['reference_pop'][0]
    target_pop = config_dict['target_pop']
    target_pop.extend([ref_pop])

    pop_file_name = 's_pop.txt'
    pop_dict = {}
    for pop in target_pop:
        pop_dict[pop] = []
    f = open(data_path + pop_file_name)
    f.readline()
    for line in f:
        line = line.replace('\n', '')
        curr_pop_data = line.split('\t')
        if curr_pop_data[1] in pop_dict:
            pop_dict[curr_pop_data[1]].append(curr_pop_data[0])
    f.close()

    reference_part = float(config_dict['reference_part'][0])
    reference_list = pop_dict[ref_pop][:int(len(pop_dict[ref_pop]) * reference_part)]
    reference_frequencies = [0, 0]

    suffix = config_dict['result_file_suffix'][0] + '_comb_' + config_dict['combinations'][0]

    result_file = open(result_path + '_'.join(target_pop) + '_' + suffix + '_result.txt', 'w')
    feature_file = open(result_path + '_'.join(target_pop) + '_' + suffix + '_features.txt', 'w')

    data_gene_list_file_name = config_dict['gene_list_name'][0]
    data_gene_file = open(data_path + data_gene_list_file_name)
    all_genes = [line.replace('\n', '') for line in data_gene_file]
    data_gene_file.close()

    L = int(config_dict['combinations'][0])

    for subset in itertools.combinations(all_genes, L):
        genes = list(subset)

        print(genes)
        result_file.write(';'.join(genes))
        result_file.write('\n')
        feature_file.write(';'.join(genes))
        feature_file.write('\n')

        main_data = []

        header = ''
        for dir_name in os.listdir(data_path):
            if dir_name.startswith('chr'):
                chr_path = data_path + dir_name + '/'
                for gene_file_name in os.listdir(chr_path):
                    if gene_file_name[:-4] in genes:
                        f = open(chr_path + gene_file_name)
                        for line in f:
                            if header == '':
                                header = line
                                main_data.append(header)
                            elif line == header:
                                continue
                            else:
                                main_data.append(line)
                        f.close()

        # Reference group

        header = main_data[0].replace('\n', '')
        header = header.split(' ')
        samples_names = header[15:]

        target_samples_ids_mt = []
        target_samples_names = []
        snp_samples_mt = {}

        number_mt_snps = 0

        for sample_name in samples_names:
            if sample_name in reference_list:
                target_samples_ids_mt.append(samples_names.index(sample_name))
                target_samples_names.append(sample_name)
                snp_samples_mt[sample_name] = []

        for item in main_data[1:]:
            number_mt_snps += 1
            item = item.replace('\n', '')
            curr_snp_data = item.split(' ')
            snp_chr = curr_snp_data[0]
            snp_data = curr_snp_data[15:]

            snp_data = list(snp_data[i] for i in target_samples_ids_mt)

            if snp_chr == 'chrMT':
                for id in range(0, len(snp_data)):
                    if snp_data[id] == '0':
                        reference_frequencies[0] += 1
                    elif snp_data[id] == '1':
                        reference_frequencies[1] += 1
        f.close()

        reference_frequencies = [freq / sum(reference_frequencies) for freq in reference_frequencies]

        # Remaining group

        header_mt = main_data[0].replace('\n', '')
        header_mt = header_mt.split(' ')
        samples_names_mt = header_mt[15:]

        target_samples_ids_mt = []
        target_samples_names_mt = []

        for sample_name in samples_names_mt:
            for pop in target_pop:
                if sample_name in pop_dict[pop]:
                    target_samples_names_mt.append(sample_name)
                    target_samples_ids_mt.append(samples_names_mt.index(sample_name))

        df_ref_mt = np.empty(shape=(len(target_samples_names_mt), number_mt_snps), dtype=float)
        names_mt = []

        line_count_mt = 0
        for item in main_data[1:]:
            if line_count_mt % 1000 == 0:
                print('Lines in mtDNA file: ', line_count_mt)

            item = item.replace('\n', '')
            curr_snp_data_mt = item.split(' ')
            snp_chr_mt = curr_snp_data_mt[0]
            snp_gene_mt = curr_snp_data_mt[1]
            snp_pos_mt = curr_snp_data_mt[2]
            snp_data_mt = curr_snp_data_mt[15:]

            snp_data_mt = list(snp_data_mt[i] for i in target_samples_ids_mt)

            if snp_chr_mt == 'chrMT':

                combination_data = []

                for id in range(0, len(snp_data_mt)):
                    if snp_data_mt[id] == '0':
                        combination_data.append(1 - reference_frequencies[0])
                    elif snp_data_mt[id] == '1':
                        combination_data.append(1 - reference_frequencies[1])

                df_ref_mt[:, line_count_mt] = combination_data

                if len(set(combination_data)) > 1:
                    name_mt = snp_gene_mt + '_' + snp_pos_mt
                    if name_mt not in names_mt:
                        names_mt.append(name_mt)

            line_count_mt += 1

        df_ref_mt = df_ref_mt[:, ~np.all(df_ref_mt[1:] == df_ref_mt[:-1], axis=0)]

        data_classes = []
        for item in target_samples_names_mt:
            for pop in target_pop:
                if item in pop_dict[pop]:
                    data_classes.append(pop)

        factor = pd.factorize(data_classes)
        y = factor[0]

        clf = RandomForestClassifier(n_estimators=10)
        output = cross_validate(clf, df_ref_mt, y, cv=5, scoring='accuracy', return_estimator=True)
        accuracy = np.mean(output['test_score'])

        print(accuracy)

        result_file.write(str(accuracy))
        result_file.write('\n')

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
        feature_file.write('Number of features: ' + str(len(features_dict.keys())))
        feature_file.write('\n')

        features_count = 0
        for key, value in features_dict.items():
            if value > 0.0 and features_count < 100:
                line = str(key) + '\t' + str(value)
                feature_file.write(line)
                feature_file.write('\n')
                features_count += 1

    result_file.close()
    feature_file.close()
