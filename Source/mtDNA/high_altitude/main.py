from Source.mtDNA.tibet.functions.file_system import get_path
from Source.mtDNA.high_altitude.functions import *

path = get_path()
tibet_data_path = path + '/Data/tibet/'
world_data_path = path + '/Data/world/'

tibet_result_path = path + '/Result/tibet/'
if not os.path.exists(tibet_result_path):
    os.makedirs(tibet_result_path)

tibet_data, tibet_subjects, tibet_classes = read_data(tibet_data_path)
world_data, world_subjects, world_classes = read_data(world_data_path)

regions = get_region_info(tibet_data_path)

current_tibet_classes = {'low': ['0-500', '501-1000', '1001-1500', '1501-2000', '2001-2500', '2501-3000'],
                         'high': ['3001-4000', '4001']}
tibet_subset, tibet_subject_classes = subset_subjects(tibet_data, tibet_classes, current_tibet_classes)
tibet_table, tibet_mutated_positions = create_classes_table(tibet_subset)

tibet_accuracy, tibet_features, tibet_accuracy_list, tibet_features_rating = \
    run_sequential_random_forest(tibet_table, tibet_subject_classes, tibet_mutated_positions, 10)

save_results(tibet_result_path, 'tibet_rf', [tibet_accuracy] + tibet_features)
save_results(tibet_result_path, 'tibet_rf_accuracy', tibet_accuracy_list)
save_results(tibet_result_path, 'tibet_rf_features', tibet_features_rating)
olo = 0
