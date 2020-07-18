from Source.mtDNA.tibet.functions.file_system import get_path
from Source.mtDNA.high_altitude.functions import *
from Source.mtDNA.high_altitude.infrastructure_functions import *

read_tibet_data = 1

path = get_path()
tibet_data_path = path + '/Data/tibet/'
world_data_path = path + '/Data/world/'

tibet_result_path = path + '/Result/tibet/'
if not os.path.exists(tibet_result_path):
    os.makedirs(tibet_result_path)

world_result_path = path + '/Result/world/'
if not os.path.exists(world_result_path):
    os.makedirs(world_result_path)

tibet_data, tibet_subjects, tibet_classes = read_data(tibet_data_path)
regions = get_region_info(tibet_data_path)

current_tibet_classes = {'low': ['0-500', '501-1000', '1001-1500', '1501-2000', '2001-2500', '2501-3000'],
                         'high': ['3001-4000', '4001']}
tibet_subset, tibet_subject_classes = subset_subjects(tibet_data, tibet_classes, current_tibet_classes)
tibet_table, tibet_mutated_positions = create_classes_table(tibet_subset)

if read_tibet_data:
    tibet_results = read_results(tibet_result_path, 'tibet_rf.txt')
    tibet_accuracy = tibet_results[0]
    tibet_features = [int(item) for item in tibet_results[1:]]
else:
    tibet_accuracy, tibet_features, tibet_accuracy_list, tibet_features_rating = \
        run_sequential_random_forest(tibet_table, tibet_subject_classes, tibet_mutated_positions, 'max')

    save_results(tibet_result_path, 'tibet_rf', [tibet_accuracy] + tibet_features)
    save_results(tibet_result_path, 'tibet_rf_accuracy', tibet_accuracy_list)
    save_results(tibet_result_path, 'tibet_rf_features', tibet_features_rating)

tibet_haplogroups = read_haplogroups(tibet_data_path, current_tibet_classes)
positions_to_remove = get_haplogroups_positions(tibet_data_path, tibet_haplogroups)
positions_to_remove_corrected = [(item - 1) for item in positions_to_remove]

tibet_filtered_features = remove_items_from_list(tibet_features, positions_to_remove_corrected)
tibet_filtered_features_corrected = [(item + 1) for item in tibet_filtered_features]

world_data, world_subjects, world_classes = read_data(world_data_path)
frequency_dict = calculate_mutation_frequency(world_data, world_classes, tibet_filtered_features_corrected)
regions_dict = calculate_regions_statistics(tibet_filtered_features_corrected, regions)

write_frequency_to_xlsx(world_result_path, 'frequency', frequency_dict)
write_regions_to_xlsx(world_result_path, 'regions', regions_dict)
