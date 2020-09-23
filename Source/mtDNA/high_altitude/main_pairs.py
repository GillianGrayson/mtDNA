from Source.mtDNA.tibet.functions.file_system import get_path
from Source.mtDNA.high_altitude.functions import *
from Source.mtDNA.high_altitude.infrastructure_functions import *

read_tibet_data = 1

path = get_path()
info_data_path = path + '/Data/alignment/info/'
tibet_data_path = path + '/Data/alignment/tibet/'
world_data_path = path + '/Data/alignment/world/'

tibet_result_path = path + '/Result/tibet_pair/'
if not os.path.exists(tibet_result_path):
    os.makedirs(tibet_result_path)

world_result_path = path + '/Result/world_pair/'
if not os.path.exists(world_result_path):
    os.makedirs(world_result_path)

tibet_data, tibet_subjects, tibet_classes = read_data(tibet_data_path)
regions = get_region_info(info_data_path)

current_tibet_classes = {
    'Asian Low Altitude': ['0-500', '501-1000', '1001-1500', '1501-2000', '2001-2500', '2501-3000', '3001-4000'],
    'Tibetan High Altitude': ['4001']}
tibet_subset, tibet_subject_classes = subset_subjects(tibet_data, tibet_classes, current_tibet_classes)

if read_tibet_data:
    tibet_results = read_results(tibet_result_path, 'tibet_rf.txt')
    tibet_accuracy = tibet_results[0]
    tibet_features = [item for item in tibet_results[1:]]
else:
    tibet_variable_positions = get_variable_positions(tibet_subset)
    tibet_table, tibet_mutated_pairs = create_classes_table_pairs(tibet_subset, tibet_variable_positions)

    tibet_accuracy, tibet_features, tibet_accuracy_list, tibet_features_rating = \
        run_sequential_random_forest(tibet_table, tibet_subject_classes, tibet_mutated_pairs, 0)

    save_results(tibet_result_path, 'tibet_rf', [tibet_accuracy] + tibet_features)
    save_results(tibet_result_path, 'tibet_rf_accuracy', tibet_accuracy_list)
    save_results(tibet_result_path, 'tibet_rf_features', tibet_features_rating)

tibet_haplogroups = read_haplogroups(info_data_path, current_tibet_classes)
positions_to_remove = get_haplogroups_positions(info_data_path, tibet_haplogroups)
positions_to_remove_corrected = [(item - 1) for item in positions_to_remove]  # haplogroups data has numeration from 1

tibet_filtered_pairs = remove_pairs_from_list(tibet_features, positions_to_remove_corrected)
tibet_filtered_pairs_items = remove_items_from_pair(tibet_features, positions_to_remove_corrected)

world_data, world_subjects, world_classes = read_data(world_data_path)
current_world_classes = {'Tibetan': ['Tibetan'], 'Andes': ['Andes'], 'Ethiopia': ['Ethiopia']}
world_subset, world_subject_classes = subset_subjects(world_data, world_classes, current_world_classes)

classes = ['Asian Low Altitude', 'Tibetan High Altitude', 'Tibetan', 'Andes', 'Ethiopia']
data = tibet_subset
data.update(world_subset)
stat_dict_pairs = create_pair_statistics(data, tibet_filtered_pairs, classes)
write_stat_to_xlsx(world_result_path, 'mutation_stat_pairs_filtered', stat_dict_pairs, 1000)

stat_dict_pairs_items = create_pair_statistics(data, tibet_filtered_pairs_items, classes)
write_stat_to_xlsx(world_result_path, 'mutation_stat_items_filtered', stat_dict_pairs_items, 1000)
