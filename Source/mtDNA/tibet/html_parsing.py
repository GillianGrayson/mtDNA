import os
import pandas as pd
from bs4 import BeautifulSoup
from Source.mtDNA.tibet.functions.file_system import get_path

path_in = get_path()
path_in += '/Data/phylotree/html/'
path_out = get_path()
path_out += '/Data/phylotree/xlsx/'

phylotrees = {'tree_name': [], 'tree': []}
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
        phylotrees['tree_name'].append(tree_name)
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
                        mutations = [elem.replace('\n', '').strip() for elem in curr_mutations]
                        tree_dict_raw['mutations'][-1].extend(mutations)
                        break
                    else:
                        tree_dict_raw['haplogroup'].append(curr_row[j])
                        tree_dict_raw['index'].append(j)
                        curr_mutations = curr_row[j + 1].split('  ')
                        mutations = [elem.replace('\n', '').strip() for elem in curr_mutations]
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
                list(dict.fromkeys([elem.strip() for j in indexes for elem in tree_dict_raw['mutations'][j]])))

        phylotrees['tree'].append(tree_dict)
        for i in range(0, len(tree_dict['mutations'])):
            mutations = tree_dict['mutations'][i]
            tree_dict['mutations'][i] = ' '.join(mutations)
        tree_df = pd.DataFrame(tree_dict)
        writer = pd.ExcelWriter(path_out + 'tree_' + tree_name + '.xlsx', engine='xlsxwriter')
        tree_df.to_excel(writer, index=False, startrow=0)
        worksheet = writer.sheets['Sheet1']
        writer.save()

phylotrees_info = {'phylotree': [], 'haplogroup': [], 'mutation_type': [], 'position': [],
                   'ancestral_base': [], 'derived_base': []}
for tree_id in range(0, len(phylotrees['tree_name'])):
    for haplogroup_id in range(0, len(phylotrees['tree'][tree_id]['haplogroup'])):
        haplogroup = phylotrees['tree'][tree_id]['haplogroup'][haplogroup_id]
        if haplogroup in phylotrees_info['haplogroup']:
            continue
        else:
            mutations = phylotrees['tree'][tree_id]['mutations'][haplogroup_id]
            mutations = mutations.split(' ')
            mutation_types = []
            positions = []
            ancestral_bases = []
            derived_bases = []
            excluded_mutations = ['reserved']
            for mutation in mutations:
                if len(mutation) == 0 or mutation in excluded_mutations:
                    continue
                if mutation.startswith('('):
                    mutation = mutation[1:-1]
                if len(mutation.split('.')) > 1:
                    mutation = mutation.split('.')
                    if mutation[1].endswith('!'):
                        excluded_mutations.append(mutation[0][1:] + '.' + mutation[1][0] + mutation[0][0])
                    else:
                        mutation_types.append('insertion')
                        positions.append(mutation[0])
                        ancestral_bases.append('')
                        try:
                            int(mutation[1][0])
                            derived_bases.append(mutation[1][1:] * int(mutation[1][0]))
                        except ValueError:
                            derived_bases.append(mutation[1])
                elif mutation.endswith('d'):
                    mutation_types.append('deletion')
                    if len(mutation[:-1].split('-')) > 1:
                        positions.append(mutation[:-1])
                        ancestral_bases.append('')
                        derived_bases.append('')
                    else:
                        positions.append(mutation[1:-1])
                        ancestral_bases.append(mutation[0])
                        derived_bases.append('')
                elif mutation.endswith('!!'):
                    excluded_mutations.append(mutation[-3] + mutation[1:-3] + mutation[0])
                elif mutation.endswith('!'):
                    excluded_mutations.append(mutation[-2] + mutation[1:-2] + mutation[0])
                    mutation_types.append('mutation')
                    positions.append(mutation[1:-2])
                    ancestral_bases.append(mutation[0])
                    derived_bases.append(mutation[-2])
                else:
                    mutation_types.append('mutation')
                    positions.append(mutation[1:-1])
                    ancestral_bases.append(mutation[0])
                    derived_bases.append(mutation[-1])
            phylotrees_info['phylotree'].extend([phylotrees['tree_name'][tree_id]] * len(positions))
            phylotrees_info['haplogroup'].extend([haplogroup] * len(positions))
            phylotrees_info['mutation_type'].extend(mutation_types)
            phylotrees_info['position'].extend(positions)
            phylotrees_info['ancestral_base'].extend(ancestral_bases)
            phylotrees_info['derived_base'].extend(derived_bases)

info_df = pd.DataFrame(phylotrees_info)
writer = pd.ExcelWriter(path_out + 'phylotrees.xlsx', engine='xlsxwriter')
info_df.to_excel(writer, index=False, startrow=0)
worksheet = writer.sheets['Sheet1']
writer.save()
