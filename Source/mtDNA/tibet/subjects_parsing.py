from Source.mtDNA.tibet.functions.file_system import get_path
import pandas as pd
import os


subject_info = {'subject': [], 'group': [], 'height': []}
data_path = get_path() + '/Data/'
for filename in os.listdir(data_path):
    if filename.endswith('.fasta'):
        f = open(data_path + filename, 'r')
        subjects_lines = [line.rstrip() for line in f][0::2]
        for subject_line in subjects_lines:
            line_list = subject_line.split(' ')
            subject_info['subject'].append(line_list[0][1:])
            group_index = line_list.index('mitochondrion,') - 1
            subject_info['group'].append(line_list[group_index])
            subject_info['height'].append(filename[:-6])
        f.close()

df = pd.DataFrame(subject_info)
writer = pd.ExcelWriter(data_path + 'subjects.xlsx', engine='xlsxwriter')
df.to_excel(writer, index=False, startrow=0)
worksheet = writer.sheets['Sheet1']
writer.save()
