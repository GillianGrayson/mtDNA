import os
import itertools
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_validate


def rf_mt_nuc(config_dict):
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
    reference_frequencies = [0, 0, 0, 0, 0, 0]

    suffix = config_dict['result_file_suffix'][0] + '_comb_' + '_'.join(config_dict['combinations'])

    result_file = open(result_path + '_'.join(target_pop) + '_' + suffix + '_result.txt', 'w')
    feature_file = open(result_path + '_'.join(target_pop) + '_' + suffix + '_features.txt', 'w')

    data_gene_list_file_name = config_dict['gene_list_name'][1]
    data_gene_file = open(data_path + data_gene_list_file_name)
    gene_list = [line.replace('\n', '') for line in data_gene_file]
    nuc_genes = []
    for dir_name in os.listdir(data_path):
        if dir_name.startswith('chr'):
            chr_path = data_path + dir_name + '/'
            for gene_file_name in os.listdir(chr_path):
                if gene_file_name[:-4] in gene_list:
                    nuc_genes.append(gene_file_name[:-4])
    data_gene_file.close()

    data_gene_list_file_name = config_dict['gene_list_name'][0]
    data_gene_file = open(data_path + data_gene_list_file_name)
    gene_list = [line.replace('\n', '') for line in data_gene_file]
    mt_genes = []
    for dir_name in os.listdir(data_path):
        if dir_name.startswith('chr'):
            chr_path = data_path + dir_name + '/'
            for gene_file_name in os.listdir(chr_path):
                if gene_file_name[:-4] in gene_list:
                    mt_genes.append(gene_file_name[:-4])
    data_gene_file.close()

    L_nuc = int(config_dict['combinations'][0])

    for nuc_subset in itertools.combinations(nuc_genes, L_nuc):
        curr_nuc_genes = list(nuc_subset)

        nuc_data = []

        header = ''
        for dir_name in os.listdir(data_path):
            if dir_name.startswith('chr'):
                chr_path = data_path + dir_name + '/'
                for gene_file_name in os.listdir(chr_path):
                    if gene_file_name[:-4] in curr_nuc_genes:
                        f = open(chr_path + gene_file_name)
                        for line in f:
                            if header == '':
                                header = line
                                nuc_data.append(header)
                            elif line == header:
                                continue
                            else:
                                nuc_data.append(line)
                        f.close()

        if len(config_dict['combinations']) == 1:
            L_mt_fin = len(mt_genes) + 1
        else:
            L_mt_fin = int(config_dict['combinations'][1]) + 1

        for L_mt in range(1, L_mt_fin):
            for mt_subset in itertools.combinations(mt_genes, L_mt):
                curr_mt_genes = list(mt_subset)

                all_genes = curr_nuc_genes + curr_mt_genes

                print(all_genes)
                result_file.write(';'.join(all_genes))
                result_file.write('\n')
                feature_file.write(';'.join(all_genes))
                feature_file.write('\n')

                mt_data = []

                header = ''
                for dir_name in os.listdir(data_path):
                    if dir_name.startswith('chr'):
                        chr_path = data_path + dir_name + '/'
                        for gene_file_name in os.listdir(chr_path):
                            if gene_file_name[:-4] in curr_mt_genes:
                                f = open(chr_path + gene_file_name)
                                for line in f:
                                    if header == '':
                                        header = line
                                        mt_data.append(header)
                                    elif line == header:
                                        continue
                                    else:
                                        mt_data.append(line)
                                f.close()

                # Reference group

                header = mt_data[0].replace('\n', '')
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

                for item in mt_data[1:]:
                    number_mt_snps += 1
                    item = item.replace('\n', '')
                    curr_snp_data = item.split(' ')
                    snp_chr = curr_snp_data[0]
                    snp_data = curr_snp_data[15:]

                    snp_data = list(snp_data[i] for i in target_samples_ids_mt)

                    if snp_chr == 'chrMT':
                        for id in range(0, len(snp_data)):
                            sample_name = samples_names[target_samples_ids_mt[id]]
                            snp_samples_mt[sample_name].append(snp_data[id])

                number_nuc_snps = 0

                header = nuc_data[0].replace('\n', '')
                header = header.split(' ')
                samples_names = header[15:]
                target_samples = []
                for sample_name in samples_names:
                    if sample_name in reference_list and sample_name in snp_samples_mt:
                        target_samples.append(samples_names.index(sample_name))

                line_count = 0
                for item in nuc_data[1:]:

                    number_nuc_snps += 1

                    if line_count % 100 == 0:
                        print('nucDNA lines: ', line_count)

                    item = item.replace('\n', '')
                    curr_snp_data = item.split(' ')
                    snp_data = curr_snp_data[15:]

                    snp_data = list(snp_data[i] for i in target_samples)

                    sample_id = 0
                    for sample in target_samples_names:
                        for id in range(0, len(snp_samples_mt[sample])):
                            if snp_samples_mt[sample][id] == '0':
                                if snp_data[sample_id] == '0|0':
                                    reference_frequencies[0] += 1
                                elif snp_data[sample_id] == '0|1' or snp_data[sample_id] == '1|0':
                                    reference_frequencies[1] += 1
                                elif snp_data[sample_id] == '1|1':
                                    reference_frequencies[2] += 1
                            elif snp_samples_mt[sample][id] == '1':
                                if snp_data[sample_id] == '0|0':
                                    reference_frequencies[3] += 1
                                elif snp_data[sample_id] == '0|1' or snp_data[sample_id] == '1|0':
                                    reference_frequencies[4] += 1
                                elif snp_data[sample_id] == '1|1':
                                    reference_frequencies[5] += 1
                        sample_id += 1
                    line_count += 1

                reference_frequencies = [freq / sum(reference_frequencies) for freq in reference_frequencies]
                snp_samples_mt = None

                # Remaining group

                header_mt = mt_data[0].replace('\n', '')
                header_mt = header_mt.split(' ')
                samples_names_mt = header_mt[15:]

                target_samples_ids_mt = []
                target_samples_names_mt = []

                for sample_name in samples_names_mt:
                    for pop in target_pop:
                        if sample_name in pop_dict[pop]:
                            target_samples_names_mt.append(sample_name)
                            target_samples_ids_mt.append(samples_names_mt.index(sample_name))

                df_ref_mt = np.empty(shape=(len(target_samples_names_mt), number_mt_snps), dtype=str)
                names_mt = []

                line_count_mt = 0
                for item in mt_data[1:]:
                    if line_count_mt % 1000 == 0:
                        print('mtDNA lines: ', line_count_mt)

                    item = item.replace('\n', '')
                    curr_snp_data_mt = item.split(' ')
                    snp_chr_mt = curr_snp_data_mt[0]
                    snp_gene_mt = curr_snp_data_mt[1]
                    snp_pos_mt = curr_snp_data_mt[2]
                    snp_data_mt = curr_snp_data_mt[15:]

                    snp_data_mt = list(snp_data_mt[i] for i in target_samples_ids_mt)

                    if snp_chr_mt == 'chrMT':

                        name_mt = snp_gene_mt + '_' + snp_pos_mt
                        if name_mt not in names_mt:
                            names_mt.append(name_mt)

                        df_ref_mt[:, line_count_mt] = snp_data_mt

                    line_count_mt += 1

                header_nuc = nuc_data[0].replace('\n', '')
                header_nuc = header_nuc.split(' ')
                samples_names_nuc = header_nuc[15:]
                target_samples_ids_nuc = []
                target_samples_names_nuc = []
                for sample_name in samples_names_nuc:
                    for pop in target_pop:
                        if sample_name in pop_dict[pop] and sample_name in target_samples_names_mt:
                            target_samples_names_nuc.append(sample_name)
                            target_samples_ids_nuc.append(samples_names_nuc.index(sample_name))

                df_ref = np.empty(shape=(len(target_samples_names_nuc), number_nuc_snps * number_mt_snps), dtype=float)
                names_combinations = []

                line_count_nuc = 0
                for item in nuc_data[1:]:

                    if line_count_nuc % 50 == 0:
                        print('Combinations: ', line_count_nuc)

                    item = item.replace('\n', '')
                    curr_snp_data_nuc = item.split(' ')
                    snp_gene_nuc = curr_snp_data_nuc[1]
                    snp_name_nuc = curr_snp_data_nuc[3]
                    snp_data_nuc = curr_snp_data_nuc[15:]

                    snp_data_nuc = list(snp_data_nuc[i] for i in target_samples_ids_nuc)

                    for mt_id in range(0, len(names_mt)):

                        combination_data = []

                        for id in range(0, len(snp_data_nuc)):
                            if df_ref_mt[id][mt_id] == '0':
                                if snp_data_nuc[id] == '0|0':
                                    combination_data.append(1 - reference_frequencies[0])
                                elif snp_data_nuc[id] == '0|1' or snp_data_nuc[id] == '1|0':
                                    combination_data.append(1 - reference_frequencies[1])
                                elif snp_data_nuc[id] == '1|1':
                                    combination_data.append(1 - reference_frequencies[2])
                            elif df_ref_mt[id][mt_id] == '1':
                                if snp_data_nuc[id] == '0|0':
                                    combination_data.append(1 - reference_frequencies[3])
                                elif snp_data_nuc[id] == '0|1' or snp_data_nuc[id] == '1|0':
                                    combination_data.append(1 - reference_frequencies[4])
                                elif snp_data_nuc[id] == '1|1':
                                    combination_data.append(1 - reference_frequencies[5])

                        df_ref[:, line_count_nuc] = combination_data

                        if len(set(combination_data)) > 1:
                            name_mt = names_mt[mt_id]
                            combination_name = name_mt + '_' + snp_gene_nuc + '_' + snp_name_nuc
                            if combination_name not in names_combinations:
                                names_combinations.append(combination_name)

                        line_count_nuc += 1

                df_ref = df_ref[:, ~np.all(df_ref[1:] == df_ref[:-1], axis=0)]

                data_classes = []
                for item in target_samples_names_nuc:
                    for pop in target_pop:
                        if item in pop_dict[pop]:
                            data_classes.append(pop)

                factor = pd.factorize(data_classes)
                y = factor[0]

                clf = RandomForestClassifier(n_estimators=10)
                output = cross_validate(clf, df_ref, y, cv=5, scoring='accuracy', return_estimator=True)
                accuracy = np.mean(output['test_score'])

                print(accuracy)

                result_file.write(str(accuracy))
                result_file.write('\n')

                features_dict = dict((key, []) for key in names_combinations)

                for idx, estimator in enumerate(output['estimator']):
                    feature_importances = pd.DataFrame(estimator.feature_importances_,
                                                       index=names_combinations,
                                                       columns=['importance']).sort_values('importance',
                                                                                           ascending=False)
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