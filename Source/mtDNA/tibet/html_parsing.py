import os
import re
import pandas as pd
from bs4 import BeautifulSoup
from mtDNA.tibet.functions.file_system import get_path

path = get_path()
path += '/Data/phylotree/'

for curr_file in os.listdir(path):
    with open(path + curr_file, "r") as f:
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
        last_group = ''
        for i in range(0, len(table_rows_list)):
            curr_row = table_rows_list[i]
            for j in range(0, len(curr_row)):
                if curr_row[j] == ' ':
                    continue
                else:
                    if curr_row[j + 1] == ' ':
                        curr_row[j - 1] = last_group
                        break
                    else:
                        last_group = curr_row[j]
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

        tree_text = soup.get_text()
        tree_text = re.sub(r'<.*?>', '', tree_text, flags=re.DOTALL)
        tree_text = re.sub(r'\n\s*\n', '\n', tree_text)
        tree_name_index = tree_text.find('subtree') + 8
        if tree_text[tree_name_index + 1] == '\n':
            tree_name = tree_text[tree_name_index]
        else:
            tree_name = tree_text[tree_name_index:tree_name_index + 1]
        tree = re.sub(r'.*mt-MRCA', 'mt-MRCA', tree_text, flags=re.DOTALL)
        tree = tree.split('var', maxsplit=1)[0]
        tree_lines = tree.split('\n')
