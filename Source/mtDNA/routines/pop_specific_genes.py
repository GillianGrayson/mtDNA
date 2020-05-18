import os
import json
import hashlib

data_path = 'E:/YandexDisk/mtDNA/Result/files/'
gene_path = 'E:/YandexDisk/mtDNA/Data/'
experiment_type = 'mt-nuc'
random_forest_type = 3
target_accuracy = 0.5
gene_files = ['mt_gene_list.txt', 'test_gene_list_cold_adaptation.txt']
reference_pops = ['GBR', 'FIN', 'TSI']
target_pops = ['GBR', 'FIN', 'TSI']
remove_genes = 0

mt_genes = {}
nuc_genes = {}
mt_nuc_genes = {}
accuracy = {}
for reference in reference_pops:
    mt_genes[reference] = {}
    nuc_genes[reference] = {}
    mt_nuc_genes[reference] = {}
    accuracy[reference] = {}
    for target in target_pops:
        if target != reference:
            mt_genes[reference][target] = []
            nuc_genes[reference][target] = []
            mt_nuc_genes[reference][target] = []
            accuracy[reference][target] = 0.0

for reference_pop in reference_pops:
    for target_pop in target_pops:
        if reference_pop != target_pop:
            genes_mt = []
            genes_names_mt = []
            genes_nuc = []
            genes_names_nuc = []
            for file_id in range(0, len(gene_files)):
                data_gene_file = open(gene_path + gene_files[file_id])
                if file_id == 0:
                    if experiment_type == 'mt':
                        for i, line in enumerate(data_gene_file):
                            genes_mt.append(i)
                            genes_names_mt.append(line.rstrip())
                        genes_nuc = []
                        genes_names_nuc = []
                    elif experiment_type == 'nuc':
                        genes_mt = []
                        genes_names_mt = []
                        for i, line in enumerate(data_gene_file):
                            genes_nuc.append(i)
                            genes_names_nuc.append(line.rstrip())
                    elif experiment_type == 'mt-nuc':
                        for i, line in enumerate(data_gene_file):
                            genes_mt.append(i)
                            genes_names_mt.append(line.rstrip())
                else:
                    for i, line in enumerate(data_gene_file):
                        genes_nuc.append(i)
                        genes_names_nuc.append(line.rstrip())
                data_gene_file.close()

            json_list = json.dumps([genes_mt, genes_nuc]).encode('utf-8')
            curr_hash = hashlib.md5(json_list).hexdigest()

            path = data_path + experiment_type + '/rf_type_' + str(random_forest_type) + \
                   '/ref_' + reference_pop + '_target_' + target_pop + '/' + curr_hash + '/'

            max_accuracy = 0.0
            max_index = 0
            accuracy_filename = str(target_accuracy) + '_accuracy.txt'
            f = open(path + accuracy_filename)
            for i, line in enumerate(f):
                curr_accuracy = float(line.rstrip())
                if curr_accuracy > max_accuracy:
                    max_accuracy = curr_accuracy
                    max_index = i
            f.close()
            accuracy[reference_pop][target_pop] = max_accuracy

            if experiment_type == 'mt' or experiment_type == 'mt-nuc':
                mt_filename = str(target_accuracy) + '_genes_mt.txt'
                f = open(path + mt_filename)
                for i, line in enumerate(f):
                    if i == max_index:
                        mt_genes[reference_pop][target_pop] = line.rstrip().split('\t')
                f.close()
            if experiment_type == 'nuc' or experiment_type == 'mt-nuc':
                nuc_filename = str(target_accuracy) + '_genes_nuc.txt'
                f = open(path + nuc_filename)
                for i, line in enumerate(f):
                    if i == max_index:
                        nuc_genes[reference_pop][target_pop] = line.rstrip().split('\t')
                f.close()

            if experiment_type == 'mt-nuc':
                for i in range(0, len(mt_genes[reference_pop][target_pop])):
                    curr_pair = mt_genes[reference_pop][target_pop][i] + ';' + nuc_genes[reference_pop][target_pop][i]
                    mt_nuc_genes[reference_pop][target_pop].append(curr_pair)

specific_mt_genes = {reference: set() for reference in reference_pops}
specific_nuc_genes = {reference: set() for reference in reference_pops}
specific_mt_nuc_genes = {reference: set() for reference in reference_pops}
for reference_pop in reference_pops:
    curr_index = reference_pops.index(reference_pop)
    if curr_index + 1 < len(target_pops):
        if experiment_type == 'mt':
            specific_mt_genes[reference_pop] = set(mt_genes[reference_pop][target_pops[curr_index + 1]])
        elif experiment_type == 'nuc':
            specific_nuc_genes[reference_pop] = set(nuc_genes[reference_pop][target_pops[curr_index + 1]])
        else:
            specific_mt_nuc_genes[reference_pop] = set(mt_nuc_genes[reference_pop][target_pops[curr_index + 1]])
    else:
        if experiment_type == 'mt':
            specific_mt_genes[reference_pop] = set(mt_genes[reference_pop][target_pops[0]])
        elif experiment_type == 'nuc':
            specific_nuc_genes[reference_pop] = set(nuc_genes[reference_pop][target_pops[0]])
        else:
            specific_mt_nuc_genes[reference_pop] = set(mt_nuc_genes[reference_pop][target_pops[0]])

for reference_pop in reference_pops:
    for target_pop in target_pops:
        if target_pop != reference_pop:
            if experiment_type == 'mt':
                specific_mt_genes[reference_pop] = specific_mt_genes[reference_pop].intersection(
                    set(mt_genes[reference_pop][target_pop]))
                specific_mt_genes[target_pop] = specific_mt_genes[target_pop].intersection(
                    set(mt_genes[reference_pop][target_pop]))
            elif experiment_type == 'nuc':
                specific_nuc_genes[reference_pop] = specific_nuc_genes[reference_pop].intersection(
                    set(nuc_genes[reference_pop][target_pop]))
                specific_nuc_genes[target_pop] = specific_nuc_genes[target_pop].intersection(
                    set(nuc_genes[reference_pop][target_pop]))
            else:
                specific_mt_nuc_genes[reference_pop] = specific_mt_nuc_genes[reference_pop].intersection(
                    set(mt_nuc_genes[reference_pop][target_pop]))
                specific_mt_nuc_genes[target_pop] = specific_mt_nuc_genes[target_pop].intersection(
                    set(mt_nuc_genes[reference_pop][target_pop]))

for reference_pop in reference_pops:
    if experiment_type == 'mt':
        specific_mt_genes[reference_pop] = list(specific_mt_genes[reference_pop])
    elif experiment_type == 'nuc':
        specific_nuc_genes[reference_pop] = list(specific_nuc_genes[reference_pop])
    else:
        specific_mt_nuc_genes[reference_pop] = list(specific_mt_nuc_genes[reference_pop])

if remove_genes:
    mt_all = []
    nuc_all = []
    mt_nuc_all = []
    for reference_pop in reference_pops:
        if experiment_type == 'mt':
            mt_all.extend(specific_mt_genes[reference_pop])
        elif experiment_type == 'nuc':
            nuc_all.extend(specific_nuc_genes[reference_pop])
        else:
            mt_nuc_all.extend(specific_mt_nuc_genes[reference_pop])

    for reference_pop in reference_pops:
        if experiment_type == 'mt':
            specific_mt_genes[reference_pop] = [item for item in specific_mt_genes[reference_pop] if
                                                mt_all.count(item) == 1]

        elif experiment_type == 'nuc':
            specific_nuc_genes[reference_pop] = [item for item in specific_nuc_genes[reference_pop] if
                                                 nuc_all.count(item) == 1]

        else:
            specific_mt_nuc_genes[reference_pop] = [item for item in specific_mt_nuc_genes[reference_pop] if
                                                    mt_nuc_all.count(item) == 1]

for reference_pop in reference_pops:
    if experiment_type == 'mt':
        for i in range(0, len(specific_mt_genes[reference_pop])):
            curr_item = int(specific_mt_genes[reference_pop][i])
            specific_mt_genes[reference_pop][i] = genes_names_mt[curr_item]

    elif experiment_type == 'nuc':
        for i in range(0, len(specific_nuc_genes[reference_pop])):
            curr_item = int(specific_nuc_genes[reference_pop][i])
            specific_nuc_genes[reference_pop][i] = genes_names_nuc[curr_item]

    else:
        for i in range(0, len(specific_mt_nuc_genes[reference_pop])):
            curr_item = specific_mt_nuc_genes[reference_pop][i].split(';')
            curr_item_mt = int(curr_item[0])
            curr_item_nuc = int(curr_item[1])
            specific_mt_nuc_genes[reference_pop][i] = genes_names_mt[curr_item_mt] + ';' + genes_names_nuc[
                curr_item_nuc]


if remove_genes:
    suffix = ''
else:
    suffix = '_all'

for reference_pop in reference_pops:
    result_path = data_path + experiment_type + '/pop_specific/ref_' + reference_pop + '/'
    if not os.path.exists(result_path):
        os.makedirs(result_path)

    if experiment_type == 'mt':
        file_name = result_path + 'mt_genes' + suffix + '.txt'
        f = open(file_name, 'w')
        for item in specific_mt_genes[reference_pop]:
            f.write(item + '\n')
        f.close()

    elif experiment_type == 'nuc':
        file_name = result_path + 'nuc_genes' + suffix + '.txt'
        f = open(file_name, 'w')
        for item in specific_nuc_genes[reference_pop]:
            f.write(item + '\n')
        f.close()

    else:
        file_name = result_path + 'mt_nuc_genes' + suffix + '.txt'
        f = open(file_name, 'w')
        for item in specific_mt_nuc_genes[reference_pop]:
            f.write(item + '\n')
        f.close()
