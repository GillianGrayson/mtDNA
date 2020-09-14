snp_pairs_dict = {}
snp_mt_dict = {}
snp_nuc_dict = {}

pop_pairs = [['GBR', 'FIN'], ['GBR', 'TSI'], ['FIN', 'TSI']]
for pop_pair in pop_pairs:
    reference_pop = pop_pair[0]
    target_pop = pop_pair[1]

    path = 'E:/YandexDisk/mtDNA/Result/files/mt-nuc/rf_type_4/ref_' + reference_pop + '_target_' + target_pop + \
           '/ceb5200f1a65cf9b742c9f4e43cbf3e4/'

    filename = '0.5_top_features_mt_nuc.txt'
    f = open(path + filename)
    line = f.readline()
    f.close()

    line_list = line.rstrip().split(';')
    line_list = line_list[:-1]
    snp_pairs = [i.rsplit(':', 1)[0] for i in line_list]
    snp_mt = [i.split('_')[0] for i in snp_pairs]
    snp_nuc = [i.split('_')[1] for i in snp_pairs]

    snp_pairs_dict[reference_pop + '_' + target_pop] = snp_pairs
    snp_mt_dict[reference_pop + '_' + target_pop] = snp_mt
    snp_nuc_dict[reference_pop + '_' + target_pop] = snp_nuc

statistics_pairs = {'common': [], 'GBR': [], 'FIN': [], 'TSI': []}
statistics_mt = {'common': [], 'GBR': [], 'FIN': [], 'TSI': []}
statistics_nuc = {'common': [], 'GBR': [], 'FIN': [], 'TSI': []}

pairs_common_intersection = list(
    set(set(snp_pairs_dict['GBR_FIN']).intersection(snp_pairs_dict['GBR_TSI'])).intersection(snp_pairs_dict['FIN_TSI']))
pairs_GBR_intersection = list(set(snp_pairs_dict['GBR_FIN']).intersection(snp_pairs_dict['GBR_TSI']))
pairs_FIN_intersection = list(set(snp_pairs_dict['GBR_FIN']).intersection(snp_pairs_dict['FIN_TSI']))
pairs_TSI_intersection = list(set(snp_pairs_dict['GBR_TSI']).intersection(snp_pairs_dict['FIN_TSI']))

mt_common_intersection = list(
    set(set(snp_mt_dict['GBR_FIN']).intersection(snp_mt_dict['GBR_TSI'])).intersection(snp_mt_dict['FIN_TSI']))
mt_GBR_intersection = list(set(snp_mt_dict['GBR_FIN']).intersection(snp_mt_dict['GBR_TSI']))
mt_FIN_intersection = list(set(snp_mt_dict['GBR_FIN']).intersection(snp_mt_dict['FIN_TSI']))
mt_TSI_intersection = list(set(snp_mt_dict['GBR_TSI']).intersection(snp_mt_dict['FIN_TSI']))

nuc_common_intersection = list(
    set(set(snp_nuc_dict['GBR_FIN']).intersection(snp_nuc_dict['GBR_TSI'])).intersection(snp_nuc_dict['FIN_TSI']))
nuc_GBR_intersection = list(set(snp_nuc_dict['GBR_FIN']).intersection(snp_nuc_dict['GBR_TSI']))
nuc_FIN_intersection = list(set(snp_nuc_dict['GBR_FIN']).intersection(snp_nuc_dict['FIN_TSI']))
nuc_TSI_intersection = list(set(snp_nuc_dict['GBR_TSI']).intersection(snp_nuc_dict['FIN_TSI']))

statistics_pairs['common'] = pairs_common_intersection
statistics_pairs['GBR'] = pairs_GBR_intersection
statistics_pairs['FIN'] = pairs_FIN_intersection
statistics_pairs['TSI'] = pairs_TSI_intersection

statistics_mt['common'] = mt_common_intersection
statistics_mt['GBR'] = mt_GBR_intersection
statistics_mt['FIN'] = mt_FIN_intersection
statistics_mt['TSI'] = mt_TSI_intersection

statistics_nuc['common'] = nuc_common_intersection
statistics_nuc['GBR'] = nuc_GBR_intersection
statistics_nuc['FIN'] = nuc_FIN_intersection
statistics_nuc['TSI'] = nuc_TSI_intersection
olo = 0
