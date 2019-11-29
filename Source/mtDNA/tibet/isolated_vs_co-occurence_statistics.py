from mtDNA.tibet.functions.tibet_functions import *
from mtDNA.tibet.functions.infrastructure_functions import *
from mtDNA.tibet.functions.plot_functions import *

use_freq = 1

data_path = 'C:/Users/User/YandexDisk/tibet/Data/'
isolated_result_path = 'C:/Users/User/YandexDisk/tibet/Result/isolated/'
co_result_path = 'C:/Users/User/YandexDisk/tibet/Result/co-occurrence/'

groups = ['0-500', '1001-1500', '1501-2000', '2001-2500', '2501-3000', '3001-4000', '4001', '501-1000']

isolated_features_data = {}
for subset in itertools.combinations(groups, 2):
    key = '_'.join(subset)
    filename = key + '_top_features.txt'
    f = open(isolated_result_path + filename, 'r')
    curr_data = [int(line.rstrip().split('\t')[0]) for line in f]
    isolated_features_data[key] = curr_data
    f.close()
filename = '8_class_top_features.txt'
f = open(isolated_result_path + filename, 'r')
curr_data = [int(line.rstrip().split('\t')[0]) for line in f]
isolated_features_data['8_class'] = curr_data
f.close()
filename = 'binary_common_features_top_accuracy.txt'
f = open(isolated_result_path + filename, 'r')
curr_data = [int(line.rstrip().split('\t')[0]) for line in f]
isolated_features_data['binary_top'] = curr_data
f.close()

co_features_data = {}
for subset in itertools.combinations(groups, 2):
    key = '_'.join(subset)
    filename = key + '_top_features.txt'
    f = open(co_result_path + filename, 'r')
    curr_data = []
    for line in f:
        curr_pair = line.rstrip().split('\t')[0]
        curr_pair_list = curr_pair.split('_')
        curr_pair_list = [int(item) for item in curr_pair_list]
        if curr_pair_list[0] in isolated_features_data[key]:
            continue
        elif curr_pair_list[1] in isolated_features_data[key]:
            continue
        else:
            curr_data.append(curr_pair)
    co_features_data[key] = curr_data
    f.close()

key = '8_class'
filename = key + '_top_features.txt'
f = open(co_result_path + filename, 'r')
curr_data = []
for line in f:
    curr_pair = line.rstrip().split('\t')[0]
    curr_pair_list = curr_pair.split('_')
    curr_pair_list = [int(item) for item in curr_pair_list]
    if curr_pair_list[0] in isolated_features_data[key]:
        continue
    elif curr_pair_list[1] in isolated_features_data[key]:
        continue
    else:
        curr_data.append(curr_pair)
co_features_data[key] = curr_data
f.close()

key = 'binary_top'
filename = 'binary_common_features_top_accuracy.txt'
f = open(co_result_path + filename, 'r')
curr_data = []
for line in f:
    curr_pair = line.rstrip().split('\t')[0]
    curr_pair_list = curr_pair.split('_')
    curr_pair_list = [int(item) for item in curr_pair_list]
    if curr_pair_list[0] in isolated_features_data[key]:
        continue
    elif curr_pair_list[1] in isolated_features_data[key]:
        continue
    else:
        curr_data.append(curr_pair)
co_features_data[key] = curr_data
f.close()

regions = get_regions(data_path)

for key in co_features_data:
    regions_intersection = create_pair_regions_stat(co_features_data[key], regions)
    file_suffix = co_result_path + key + '_ex_regions.txt'
    save_dict(regions_intersection, file_suffix)
    plot_hist(regions_intersection, key + '_ex_regions', co_result_path)
