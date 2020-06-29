from Source.mtDNA.tibet.functions.tibet_functions import *
from Source.mtDNA.tibet.functions.infrastructure_functions import *
import pandas as pd
from Source.mtDNA.tibet.functions.file_system import get_path


path = get_path()
data_path = path + '/Data/'
result_path = path + '/Result/haplogroups/'

if not os.path.exists(result_path):
    os.makedirs(result_path)

haplogroups = get_haplogroups(data_path)

haplogroups_stat = get_haplogroups_statistics(haplogroups)

haplogroups_list = []
altitude_groups_list = []
haplogroups_counts = []
for group in haplogroups_stat:
    for key in haplogroups_stat[group]:
        haplogroups_list.append(key)
        altitude_groups_list.append(group)
        haplogroups_counts.append(haplogroups_stat[group][key])

common_haplogroups_list = []
common_altitude_groups_list = []
common_haplogroups_counts = []
for item in haplogroups_list:
    if haplogroups_list.count(item) > 1:
        indexes = [i for i, elem in enumerate(haplogroups_list) if elem == item]
        common_haplogroups_list.extend([haplogroups_list[index] for index in indexes])
        common_altitude_groups_list.extend([altitude_groups_list[index] for index in indexes])
        common_haplogroups_counts.extend([haplogroups_counts[index] for index in indexes])

haplogroups_data = {'Haplogroups': haplogroups_list,
                    'Altitude groups': altitude_groups_list,
                    'Count': haplogroups_counts}

common_haplogroups_data = {'Haplogroups': common_haplogroups_list,
                           'Altitude groups': common_altitude_groups_list,
                           'Count': common_haplogroups_counts}

haplogroups_df = pd.DataFrame(haplogroups_data)
haplogroups_df.to_excel(result_path + 'haplogroups_stat.xlsx', index=False)

common_haplogroups_df = pd.DataFrame(common_haplogroups_data)
common_haplogroups_df.to_excel(result_path + 'common_haplogroups_stat.xlsx', index=False)
