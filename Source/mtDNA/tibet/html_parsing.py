import os
from bs4 import BeautifulSoup
from mtDNA.tibet.functions.file_system import get_path

path = get_path()
path += '/Data/phylotree/'

for curr_file in os.listdir(path):
    with open(path + curr_file, "r") as f:
        contents = f.read()
        soup = BeautifulSoup(contents, 'lxml')