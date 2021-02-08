from Source.mtDNA.tibet.functions.file_system import get_path
from Source.mtDNA.high_altitude.functions import *
from Source.mtDNA.high_altitude.infrastructure_functions import *

path = get_path()
info_data_path = path + '/Data/world/rcrs/info/'
world_data_path = path + '/Data/world/rcrs/wo_hg/'

world_result_path = path + '/Result/world/rcrs/wo_hg/'
if not os.path.exists(world_result_path):
    os.makedirs(world_result_path)

world_data, world_subjects, world_classes = read_data(world_data_path)
regions = get_region_info(info_data_path)

current_world_classes = {'Tibetan': ['Tibetan'], 'Andes': ['Andes']}
world_subset, world_subject_classes = subset_subjects(world_data, world_classes, current_world_classes)
world_table, world_mutated_positions = create_classes_table(world_subset)

world_accuracy, world_features, world_accuracy_list, world_features_rating = \
    run_sequential_random_forest(world_table, world_subject_classes, world_mutated_positions, 'max')

save_results(world_result_path, 'world_rf_TA', [world_accuracy] + world_features)
save_results(world_result_path, 'world_rf_accuracy_TA', world_accuracy_list)
save_results(world_result_path, 'world_rf_features_TA', world_features_rating)

current_world_classes = {'Andes': ['Andes'], 'Ethiopia': ['Ethiopia']}
world_subset, world_subject_classes = subset_subjects(world_data, world_classes, current_world_classes)
world_table, world_mutated_positions = create_classes_table(world_subset)

world_accuracy, world_features, world_accuracy_list, world_features_rating = \
    run_sequential_random_forest(world_table, world_subject_classes, world_mutated_positions, 'max')

save_results(world_result_path, 'world_rf_AE', [world_accuracy] + world_features)
save_results(world_result_path, 'world_rf_accuracy_AE', world_accuracy_list)
save_results(world_result_path, 'world_rf_features_AE', world_features_rating)

current_world_classes = {'Tibetan': ['Tibetan'], 'Ethiopia': ['Ethiopia']}
world_subset, world_subject_classes = subset_subjects(world_data, world_classes, current_world_classes)
world_table, world_mutated_positions = create_classes_table(world_subset)

world_accuracy, world_features, world_accuracy_list, world_features_rating = \
    run_sequential_random_forest(world_table, world_subject_classes, world_mutated_positions, 'max')

save_results(world_result_path, 'world_rf_TE', [world_accuracy] + world_features)
save_results(world_result_path, 'world_rf_accuracy_TE', world_accuracy_list)
save_results(world_result_path, 'world_rf_features_TE', world_features_rating)

current_world_classes = {'Tibetan': ['Tibetan'], 'Andes': ['Andes'], 'Ethiopia': ['Ethiopia']}
world_subset, world_subject_classes = subset_subjects(world_data, world_classes, current_world_classes)
world_table, world_mutated_positions = create_classes_table(world_subset)

world_accuracy, world_features, world_accuracy_list, world_features_rating = \
    run_sequential_random_forest(world_table, world_subject_classes, world_mutated_positions, 'max')

save_results(world_result_path, 'world_rf_TAE', [world_accuracy] + world_features)
save_results(world_result_path, 'world_rf_accuracy_TAE', world_accuracy_list)
save_results(world_result_path, 'world_rf_features_TAE', world_features_rating)
