from mtDNA.tibet.functions.tibet_functions import *
from mtDNA.tibet.functions.infrastructure_functions import *
from mtDNA.tibet.functions.plot_functions import *


use_freq = 1

data_path = 'C:/Users/User/YandexDisk/tibet/Data/'
data_seq_path = 'C:/Users/User/YandexDisk/tibet/Result/co-occurrence/'
result_path = 'C:/Users/User/YandexDisk/tibet/Result/co-occurrence/seq/'

if not os.path.exists(result_path):
    os.makedirs(result_path)

raw_data, data_classes = get_data(data_path)

regions = get_regions(data_path)

variable_positions = get_variable_positions(raw_data)

if use_freq:
    df_main, position_pairs = create_co_df_freq(raw_data, variable_positions)
else:
    df_main, position_pairs = create_co_df(raw_data, variable_positions)

classes = []
for i in range(0, len(raw_data)):
    classes += [data_classes[i], ] * len(raw_data[i])

features_filename = data_seq_path + '8_class_top_features.txt'
features_dict = get_features_dict(features_filename)

top_accuracy, top_features = run_sequential_random_forest_preset(df_main, classes, position_pairs, features_dict)
print('8-class Classification Accuracy: ' + str(top_accuracy))

result_file_name = 'classification.txt'
f = open(result_path + result_file_name, 'w')
f.write('8-class Classification Accuracy: ' + str(top_accuracy) + '\n')

top_regions = create_pair_regions_stat(top_features, regions)

file_suffix = result_path + '8_class_top_features.txt'
save_list(top_features, file_suffix)

file_suffix = result_path + '8_class_regions.txt'
save_dict(top_regions, file_suffix)

plot_hist(top_regions, '8_class_regions', result_path)

features_intersection = position_pairs
features_intersection_acc = position_pairs

for subset in itertools.combinations(raw_data, 2):
    test_data = subset
    variable_positions = get_variable_positions(test_data)
    if use_freq:
        df_main, position_pairs = create_co_df_freq(test_data, variable_positions)
    else:
        df_main, position_pairs = create_co_df(test_data, variable_positions)

    classes = [data_classes[raw_data.index(test_data[0])], ] * len(test_data[0]) + \
              [data_classes[raw_data.index(test_data[1])], ] * len(test_data[1])

    features_filename = data_seq_path + data_classes[raw_data.index(test_data[0])] + '_' + \
                        data_classes[raw_data.index(test_data[1])] + '_top_features.txt'
    features_dict = get_features_dict(features_filename)

    top_accuracy, top_features = run_sequential_random_forest_preset(df_main, classes, position_pairs, features_dict)
    print(data_classes[raw_data.index(test_data[0])] + ' vs ' + data_classes[raw_data.index(test_data[1])] +
          ' Binary Classification Accuracy: ' + str(top_accuracy))
    f.write(data_classes[raw_data.index(test_data[0])] + ' vs ' + data_classes[raw_data.index(test_data[1])] +
            ' Binary Classification Accuracy: ' + str(top_accuracy) + '\n')

    top_regions = create_pair_regions_stat(top_features, regions)

    features_intersection = list(set(features_intersection).intersection(set(top_features)))
    if top_accuracy > 0.8:
        features_intersection_acc = list(set(features_intersection_acc).intersection(set(top_features)))

    file_suffix = result_path + data_classes[raw_data.index(test_data[0])] + '_' + \
                  data_classes[raw_data.index(test_data[1])] + '_top_features.txt'
    save_list(top_features, file_suffix)

    file_suffix = result_path + data_classes[raw_data.index(test_data[0])] + '_' + \
                  data_classes[raw_data.index(test_data[1])] + '_regions.txt'
    save_dict(top_regions, file_suffix)

    plot_hist(top_regions, data_classes[raw_data.index(test_data[0])] + '_' +
              data_classes[raw_data.index(test_data[1])] + '_regions', result_path)
f.close()

file_suffix = result_path + 'binary_common_features.txt'
save_list(features_intersection, file_suffix)

file_suffix = result_path + 'binary_common_features_top_accuracy.txt'
save_list(features_intersection_acc, file_suffix)

regions_intersection = create_pair_regions_stat(features_intersection, regions)
file_suffix = result_path + 'binary_common_regions.txt'
save_dict(regions_intersection, file_suffix)
plot_hist(regions_intersection, 'binary_common_regions', result_path)

regions_intersection_acc = create_pair_regions_stat(features_intersection_acc, regions)
file_suffix = result_path + 'binary_common_regions_top_accuracy.txt'
save_dict(regions_intersection_acc, file_suffix)
plot_hist(regions_intersection_acc, 'binary_common_regions_top_accuracy', result_path)
