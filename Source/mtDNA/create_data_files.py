import os
import numpy as np

data_path = '../../Data/genes/'

for file_name in os.listdir(data_path):
    with open(data_path + file_name) as f:
        for i, l in enumerate(f):
            pass
    num_lines = i
    f = open(data_path + file_name)
    header = f.readline()
    header.replace('\n', '')
    subjects = header.split(' ')[15:]
    data = np.empty(shape=(num_lines, len(subjects)), dtype=int)
    line_id = 0
    for line in f:
        line = line.replace('\n', '')
        curr_data = line.split(' ')
        chr = curr_data[0]
        gene = curr_data[1]
        pos = curr_data[2]
        name = curr_data[3]
        subject_id = 0
        for item in curr_data[15:]:
            if chr == 'chrMT':
                if item == '0':
                    data[line_id, subject_id] = 0
                elif item == '1':
                    data[line_id, subject_id] = 1
            else:
                if item == '0|0':
                    data[line_id, subject_id] = 0
                elif item == '0|1':
                    data[line_id, subject_id] = 1
                elif item == '1|0':
                    data[line_id, subject_id] = 2
                elif item == '1|1':
                    data[line_id, subject_id] = 3
            subject_id += 1
        line_id += 1
    f.close()



