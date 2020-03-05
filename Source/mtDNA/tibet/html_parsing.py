import os
import re
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

        tree_dict = {'haplogroup': [], 'index': [], 'mutations': []}
        for i in range(0, len(table_rows_list)):
            curr_row = table_rows_list[i][:-2]
            for j in range(0, len(curr_row)):
                if curr_row[j] == ' ':
                    continue
                else:
                    if curr_row[j + 1] == ' ':
                        curr_mutations = curr_row[j].split('  ')
                        mutations = [elem.replace('\n', '') for elem in curr_mutations]
                        tree_dict['mutations'][-1].extend(mutations)
                        break
                    else:
                        tree_dict['haplogroup'].append(curr_row[j])
                        tree_dict['index'].append(j)
                        curr_mutations = curr_row[j + 1].split('  ')
                        mutations = [elem.replace('\n', '') for elem in curr_mutations]
                        tree_dict['mutations'].append(mutations)
                        if j > 0:
                            tree_dict['mutations'][-1].extend(tree_dict['mutations'][j - 1])
                        break




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
