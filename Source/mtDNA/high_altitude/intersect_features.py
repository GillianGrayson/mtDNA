import os
from Source.mtDNA.tibet.functions.file_system import get_path

data = {}
path = get_path()
path += '/Result/align/01/'

for file in os.listdir(path):
    if file.endswith('txt'):
        f = open(path + file)
        f.readline()
        data[file] = []
        for line in f:
            data[file].append(line.rstrip())
        f.close()

common_data = data[list(data.keys())[0]]
for key in data:
    common_data = set(common_data).intersection(set(data[key]))

if len(common_data) > 0:
    f = open(path + 'intersection.txt', 'w')
    for item in common_data:
        f.write(item + '\n')
    f.close()
