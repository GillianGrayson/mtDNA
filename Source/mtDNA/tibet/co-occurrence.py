from mtDNA.tibet.tibet_functions import *
from mtDNA.tibet.infrastructure_functions import *
from mtDNA.tibet.plot_functions import *

use_freq = 0

data_path = 'C:/Users/User/YandexDisk/tibet/Data/'
result_path = 'C:/Users/User/YandexDisk/tibet/Result/co-occurrence/'

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

accuracy, features_dict = run_random_forest(df_main, classes, position_pairs)
print('8-class Classification Accuracy: ' + str(accuracy))

result_file_name = 'classification.txt'
f = open(result_path + result_file_name, 'w')
f.write('8-class Classification Accuracy: ' + str(accuracy) + '\n')

top_features = create_features_top(features_dict)
top_regions = create_pair_regions_stat(list(top_features.keys()), regions)

file_suffix = result_path + '8_class_top_features.txt'
save_dict(top_features, file_suffix)

file_suffix = result_path + '8_class_regions.txt'
save_dict(top_regions, file_suffix)

plot_hist(top_regions, '8_class_regions', result_path)

features_intersection = position_pairs
features_intersection_acc = position_pairs

for subset in itertools.combinations(raw_data, 2):
    test_data = subset
    if use_freq:
        df_main, position_pairs = create_co_df_freq(test_data, variable_positions)
    else:
        df_main, position_pairs = create_co_df(test_data, variable_positions)

    classes = [data_classes[raw_data.index(test_data[0])], ] * len(test_data[0]) + \
              [data_classes[raw_data.index(test_data[1])], ] * len(test_data[1])

    accuracy, features_dict = run_random_forest(df_main, classes, position_pairs)
    print(data_classes[raw_data.index(test_data[0])] + ' vs ' + data_classes[raw_data.index(test_data[1])] +
          ' Binary Classification Accuracy: ' + str(accuracy))
    f.write(data_classes[raw_data.index(test_data[0])] + ' vs ' + data_classes[raw_data.index(test_data[1])] +
            ' Binary Classification Accuracy: ' + str(accuracy) + '\n')

    top_features = create_features_top(features_dict)
    top_regions = create_pair_regions_stat(list(top_features.keys()), regions)

    features_intersection = list(set(features_intersection).intersection(list(set(top_features.keys()))))
    if accuracy > 0.8:
        features_intersection_acc = list(set(features_intersection_acc).intersection(list(set(top_features.keys()))))

    file_suffix = result_path + data_classes[raw_data.index(test_data[0])] + '_' + \
                  data_classes[raw_data.index(test_data[1])] + '_top_features.txt'
    save_dict(top_features, file_suffix)

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
