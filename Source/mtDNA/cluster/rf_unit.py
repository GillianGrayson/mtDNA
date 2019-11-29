from rf_mt import rf_type_0_mt, rf_type_1_mt, rf_type_2_mt, rf_type_3_mt
from rf_nuc import rf_type_0_nuc, rf_type_1_nuc, rf_type_2_nuc, rf_type_3_nuc
from rf_mt_nuc import rf_type_0_mt_nuc, rf_type_1_mt_nuc, rf_type_2_mt_nuc, rf_type_3_mt_nuc


def unit_task(config, results):
    if config.params_dict['experiment_type'] == 'mt':
        task_mt(config, results)
    if config.params_dict['experiment_type'] == 'nuc':
        task_nuc(config, results)
    if config.params_dict['experiment_type'] == 'mt-nuc':
        task_mt_nuc(config, results)


def task_mt(config, results):
    if int(config.params_dict['random_forest_type']) == 0:
        rf_type_0_mt(config, results)
    elif int(config.params_dict['random_forest_type']) == 1:
        rf_type_1_mt(config, results)
    elif int(config.params_dict['random_forest_type']) == 2:
        rf_type_2_mt(config, results)
    elif int(config.params_dict['random_forest_type']) == 3:
        rf_type_3_mt(config, results)


def task_nuc(config, results):
    if int(config.params_dict['random_forest_type']) == 0:
        rf_type_0_nuc(config, results)
    elif int(config.params_dict['random_forest_type']) == 1:
        rf_type_1_nuc(config, results)
    elif int(config.params_dict['random_forest_type']) == 2:
        rf_type_2_nuc(config, results)
    elif int(config.params_dict['random_forest_type']) == 3:
        rf_type_3_nuc(config, results)


def task_mt_nuc(config, results):
    if int(config.params_dict['random_forest_type']) == 0:
        rf_type_0_mt_nuc(config, results)
    elif int(config.params_dict['random_forest_type']) == 1:
        rf_type_1_mt_nuc(config, results)
    elif int(config.params_dict['random_forest_type']) == 2:
        rf_type_2_mt_nuc(config, results)
    elif int(config.params_dict['random_forest_type']) == 3:
        rf_type_3_mt_nuc(config, results)
