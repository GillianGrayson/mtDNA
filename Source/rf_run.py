from rf_mt import rf_mt
from rf_nuc import rf_nuc
from rf_mt_nuc import rf_mt_nuc

config_file_name = 'config.txt'
f = open(config_file_name)
config_dict = {}
for line in f:
    line = line.replace('\n', '')
    item = line.split(': ')
    key = item[0]
    value = item[1].split(', ')
    config_dict[key] = value
f.close()

if config_dict['experiment'][0] == 'mt':
    rf_mt(config_dict)
if config_dict['experiment'][0] == 'nuc':
    rf_nuc(config_dict)
if config_dict['experiment'][0] == 'mt-nuc':
    rf_mt_nuc(config_dict)