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
