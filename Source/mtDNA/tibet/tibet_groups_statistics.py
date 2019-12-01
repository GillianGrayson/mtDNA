from mtDNA.tibet.functions.tibet_functions import *
from mtDNA.tibet.functions.infrastructure_functions import *
from mtDNA.tibet.functions.plot_functions import *
from mtDNA.tibet.functions.file_system import get_path


use_freq = 1

path = get_path()
data_path = path + '/Data/'
result_path = path + '/Result/isolated/'

if not os.path.exists(result_path):
    os.makedirs(result_path)

groups = ['0-500', '1001-1500', '1501-2000', '2001-2500', '2501-3000', '3001-4000', '4001', '501-1000']

features_data = {key: [] for key in groups}
for subset in itertools.combinations(groups, 2):
    filename = '_'.join(subset) + '_top_features.txt'
    f = open(result_path + filename, 'r')
    curr_data = [int(line.rstrip().split('\t')[0]) for line in f]
    features_data[subset[0]].append(curr_data)
    features_data[subset[1]].append(curr_data)
    f.close()

regions = get_regions(data_path)

top_regions = {key: [] for key in groups}
for key in groups:
    curr_features = features_data[key][0]
    for i in range(1, len(features_data[key])):
        curr_features = list(set(curr_features).intersection(features_data[key][i]))
    top_regions[key] = create_regions_stat(curr_features, regions)

    file_suffix = result_path + key + '_regions.txt'
    save_dict(top_regions[key], file_suffix)

    plot_hist(top_regions[key], key + '_regions', result_path)
