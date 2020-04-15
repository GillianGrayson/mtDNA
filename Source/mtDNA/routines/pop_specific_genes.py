import json
import hashlib

data_path = 'E:/YandexDisk/mtDNA/Result/files/'
gene_path = 'E:/YandexDisk/mtDNA/Data/'
experiment_type = 'mt-nuc'
random_forest_type = 3
target_accuracy = 0.55
gene_files = ['mt_gene_list.txt', 'test_gene_list_cold_adaptation.txt']
reference_pops = ['GBR', 'FIN', 'TSI', 'IBS']
target_pops = ['GBR', 'FIN', 'TSI', 'IBS']

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
            genes_nuc = []
            for file_id in range(0, len(gene_files)):
                data_gene_file = open(gene_path + gene_files[file_id])
                if file_id == 0:
                    if experiment_type == 'mt':
                        for i, line in enumerate(data_gene_file):
                            genes_mt.append(i)
                        genes_nuc = []
                    elif experiment_type == 'nuc':
                        genes_mt = []
                        for i, line in enumerate(data_gene_file):
                            genes_nuc.append(i)
                    elif experiment_type == 'mt-nuc':
                        for i, line in enumerate(data_gene_file):
                            genes_mt.append(i)
                else:
                    for i, line in enumerate(data_gene_file):
                        genes_nuc.append(i)
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

            for i in range(0, len(mt_genes[reference_pop][target_pop])):
                curr_pair = mt_genes[reference_pop][target_pop][i] + ';' + nuc_genes[reference_pop][target_pop][i]
                mt_nuc_genes[reference_pop][target_pop].append(curr_pair)

specific_mt_genes = {reference: [] for reference in reference_pops}
specific_nuc_genes = {reference: [] for reference in reference_pops}
specific_mt_nuc_genes = {reference: [] for reference in reference_pops}
for reference_pop in reference_pops:
    for target_pop in target_pops:
        if target_pop != reference_pop:
            for item in mt_genes[reference_pop][target_pop]:
                if item in specific_mt_genes[reference_pop]:
                    specific_mt_genes[reference_pop].pop(item)
                else:
                    specific_mt_genes[reference_pop].append(item)
            for item in nuc_genes[reference_pop][target_pop]:
                if item in specific_nuc_genes[reference_pop]:
                    specific_nuc_genes[reference_pop].pop(item)
                else:
                    specific_nuc_genes[reference_pop].append(item)
            for item in mt_nuc_genes[reference_pop][target_pop]:
                if item in specific_mt_nuc_genes[reference_pop]:
                    specific_mt_nuc_genes[reference_pop].pop(item)
                else:
                    specific_mt_nuc_genes[reference_pop].append(item)
olo = 8
