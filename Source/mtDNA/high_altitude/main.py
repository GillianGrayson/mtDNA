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
tibet_data = [[tibet_data[group_id][subject_id][1:] for subject_id in range(0, len(tibet_data[group_id]))] for group_id
              in range(0, len(tibet_data))]
regions = get_region_info(tibet_data_path)

current_tibet_classes = {
    'Tibetan Low Altitude': ['0-500', '501-1000', '1001-1500', '1501-2000', '2001-2500', '2501-3000'],
    'Tibetan High Altitude': ['3001-4000', '4001']}
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
positions_to_remove_corrected = [(item - 1) for item in positions_to_remove]  # haplogroups data has numeration from 1

tibet_filtered_features = remove_items_from_list(tibet_features, positions_to_remove_corrected)

world_data, world_subjects, world_classes = read_data(world_data_path)
current_world_classes = {'Tibetan': ['Tibetan'], 'Andes': ['Andes'], 'Ethiopia': ['Ethiopia']}
world_subset, world_subject_classes = subset_subjects(world_data, world_classes, current_world_classes)

frequency_dict = calculate_mutation_frequency(world_data, world_classes, tibet_filtered_features)
frequency_dict_non_zero_1, frequency_dict_non_zero_2, frequency_dict_non_zero_3 = filter_frequency_dict(frequency_dict)
tibet_filtered_features_corrected = [(tibet_filtered_features[i] + 1) for i in range(0, len(tibet_filtered_features))]
frequency_non_zero_1_corrected = [(int(list(frequency_dict_non_zero_1.keys())[i]) + 1) for i in
                                  range(0, len(list(frequency_dict_non_zero_1.keys())))]
frequency_non_zero_2_corrected = [(int(list(frequency_dict_non_zero_2.keys())[i]) + 1) for i in
                                  range(0, len(list(frequency_dict_non_zero_2.keys())))]
frequency_non_zero_3_corrected = [(int(list(frequency_dict_non_zero_3.keys())[i]) + 1) for i in
                                  range(0, len(list(frequency_dict_non_zero_3.keys())))]
regions_dict = calculate_regions_statistics(tibet_filtered_features_corrected, regions)
regions_dict_non_zero_1 = calculate_regions_statistics(frequency_non_zero_1_corrected, regions)
regions_dict_non_zero_2 = calculate_regions_statistics(frequency_non_zero_2_corrected, regions)
regions_dict_non_zero_3 = calculate_regions_statistics(frequency_non_zero_3_corrected, regions)

write_frequency_to_xlsx(world_result_path, 'frequency', frequency_dict)
write_frequency_to_xlsx(world_result_path, 'frequency_non_zero_1', frequency_dict_non_zero_1)
write_frequency_to_xlsx(world_result_path, 'frequency_non_zero_2', frequency_dict_non_zero_2)
write_frequency_to_xlsx(world_result_path, 'frequency_non_zero_3', frequency_dict_non_zero_3)
write_regions_to_xlsx(world_result_path, 'regions', regions_dict)
write_regions_to_xlsx(world_result_path, 'regions_non_zero_1', regions_dict_non_zero_1)
write_regions_to_xlsx(world_result_path, 'regions_non_zero_2', regions_dict_non_zero_2)
write_regions_to_xlsx(world_result_path, 'regions_non_zero_3', regions_dict_non_zero_3)

classes = ['Tibetan Low Altitude', 'Tibetan High Altitude', 'Tibetan', 'Andes', 'Ethiopia']
data = tibet_subset
data.update(world_subset)
stat_dict = create_mutation_statistics(data, tibet_filtered_features, classes)
write_stat_to_xlsx(world_result_path, 'mutation_stat', stat_dict)
