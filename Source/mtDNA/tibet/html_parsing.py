import os
import pandas as pd
from bs4 import BeautifulSoup
from mtDNA.tibet.functions.file_system import get_path

path_in = get_path()
path_in += '/Data/phylotree/html/'
path_out = get_path()
path_out += '/Data/phylotree/xlsx/'

for curr_file in os.listdir(path_in):
    with open(path_in + curr_file, "r") as f:
        contents = f.read()
        soup = BeautifulSoup(contents, 'lxml')

        table = soup.find('table')
        table_rows = table.find_all('tr')
        table_rows_list = []
        for tr in table_rows:
            td = tr.find_all('td')
            row = [i.text for i in td]
            table_rows_list.append(row)

        table_rows_list = [row for row in table_rows_list if len(set(row)) > 1]
        tree_name_index = table_rows_list[0][1].find('subtree') + len('subtree ')
        tree_name = table_rows_list[0][1][tree_name_index:]
        table_rows_list = table_rows_list[16:]
        table_rows_list = [[elem.replace('\xa0', ' ') for elem in table_row] for table_row in table_rows_list]
        df = pd.DataFrame(table_rows_list)

        table_rows_list_mod = []
        groups = []
        group_indexes = []
        for i in range(0, len(table_rows_list)):
            curr_row = table_rows_list[i]
            for j in range(0, len(curr_row)):
                if curr_row[j] == ' ':
                    continue
                else:
                    if curr_row[j + 1] == ' ':
                        if j > group_indexes[-1]:
                            curr_row[j - 1] = groups[-1]
                        else:
                            nearest_index = len(group_indexes) - 1 - group_indexes[::-1].index(j - 1)
                            curr_row[j - 1] = groups[nearest_index]
                        groups.append(curr_row[j - 1])
                        group_indexes.append(j - 1)
                        break
                    else:
                        groups.append(curr_row[j])
                        group_indexes.append(j)
                        break
            table_rows_list_mod.append(curr_row)
        df_mod = pd.DataFrame(table_rows_list_mod)

        tree_dict_raw = {'haplogroup': [], 'index': [], 'mutations': []}
        for i in range(0, len(table_rows_list_mod)):
            curr_row = table_rows_list_mod[i][:-2]
            for j in range(0, len(curr_row)):
                if curr_row[j] == ' ':
                    continue
                else:
                    if curr_row[j + 1] == ' ':
                        curr_mutations = curr_row[j].split('  ')
                        mutations = [elem.replace('\n', '') for elem in curr_mutations]
                        tree_dict_raw['mutations'][-1].extend(mutations)
                        break
                    else:
                        tree_dict_raw['haplogroup'].append(curr_row[j])
                        tree_dict_raw['index'].append(j)
                        curr_mutations = curr_row[j + 1].split('  ')
                        mutations = [elem.replace('\n', '') for elem in curr_mutations]
                        tree_dict_raw['mutations'].append(mutations)
                        if j > 0:
                            nearest_index = len(tree_dict_raw['index']) - 1 - tree_dict_raw['index'][::-1].index(j - 1)
                            tree_dict_raw['mutations'][-1].extend(tree_dict_raw['mutations'][nearest_index])
                        break

        haplogroups = list(dict.fromkeys(tree_dict_raw['haplogroup']))
        tree_dict = {'haplogroup': [], 'mutations': []}
        for i in range(0, len(haplogroups)):
            haplogroup = haplogroups[i]
            indexes = [i for i, v in enumerate(tree_dict_raw['haplogroup']) if v == haplogroup]
            tree_dict['haplogroup'].append(haplogroup)
            tree_dict['mutations'].append(
                list(dict.fromkeys([elem for j in indexes for elem in tree_dict_raw['mutations'][j]])))

        for i in range(0, len(tree_dict['mutations'])):
            mutations = tree_dict['mutations'][i]
            tree_dict['mutations'][i] = ' '.join(mutations)
        tree_df = pd.DataFrame(tree_dict)
        writer = pd.ExcelWriter(path_out + 'tree_' + tree_name + '.xlsx', engine='xlsxwriter')
        tree_df.to_excel(writer, index=False, startrow=0)
        worksheet = writer.sheets['Sheet1']
        writer.save()
