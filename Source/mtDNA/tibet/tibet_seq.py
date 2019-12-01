from mtDNA.tibet.functions.tibet_functions import *
from mtDNA.tibet.functions.infrastructure_functions import *
from mtDNA.tibet.functions.plot_functions import *
from mtDNA.tibet.functions.file_system import get_path


use_freq = 1

path = get_path()
data_path = path + '/Data/'
result_path = path + '/isolated/seq/'

if not os.path.exists(result_path):
    os.makedirs(result_path)

raw_data, data_classes = get_data(data_path)

regions = get_regions(data_path)

if use_freq:
    df_main, positions = create_df_freq(raw_data)
else:
    df_main, positions = create_df(raw_data)

classes = []
for i in range(0, len(raw_data)):
    classes += [data_classes[i], ] * len(raw_data[i])

top_accuracy, top_features = run_sequential_random_forest(df_main, classes, positions)
print('8-class Classification Accuracy: ' + str(top_accuracy))

result_file_name = 'classification.txt'
f = open(result_path + result_file_name, 'w')
f.write('8-class Classification Accuracy: ' + str(top_accuracy) + '\n')

top_regions = create_regions_stat(top_features, regions)

file_suffix = result_path + '8_class_top_features.txt'
save_list(top_features, file_suffix)

file_suffix = result_path + '8_class_regions.txt'
save_dict(top_regions, file_suffix)

plot_hist(top_regions, '8_class_regions', result_path)

features_intersection = positions
features_intersection_acc = positions

for subset in itertools.combinations(raw_data, 2):
    test_data = subset
    if use_freq:
        df_main, positions = create_df_freq(test_data)
    else:
        df_main, positions = create_df(test_data)

    classes = [data_classes[raw_data.index(test_data[0])], ] * len(test_data[0]) + \
              [data_classes[raw_data.index(test_data[1])], ] * len(test_data[1])

    top_accuracy, top_features = run_sequential_random_forest(df_main, classes, positions)
    print(data_classes[raw_data.index(test_data[0])] + ' vs ' + data_classes[raw_data.index(test_data[1])] +
          ' Binary Classification Accuracy: ' + str(top_accuracy))
    f.write(data_classes[raw_data.index(test_data[0])] + ' vs ' + data_classes[raw_data.index(test_data[1])] +
            ' Binary Classification Accuracy: ' + str(top_accuracy) + '\n')

    top_regions = create_regions_stat(top_features, regions)

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

regions_intersection = create_regions_stat(features_intersection, regions)
file_suffix = result_path + 'binary_common_regions.txt'
save_dict(regions_intersection, file_suffix)
plot_hist(regions_intersection, 'binary_common_regions', result_path)

regions_intersection_acc = create_regions_stat(features_intersection_acc, regions)
file_suffix = result_path + 'binary_common_regions_top_accuracy.txt'
save_dict(regions_intersection_acc, file_suffix)
plot_hist(regions_intersection_acc, 'binary_common_regions_top_accuracy', result_path)
