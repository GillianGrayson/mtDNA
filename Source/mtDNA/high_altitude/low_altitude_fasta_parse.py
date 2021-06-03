from Source.mtDNA.tibet.functions.file_system import get_path
import os

path = get_path()
data_path = path + '/Data/low_altitude/'

subjects_file = 'low.txt'
subjects = []
f = open(data_path + subjects_file, 'r')
for line in f:
    curr_subject = line.rstrip().split('\t')[0]
    subjects.append(curr_subject)
f.close()

data = []
for curr_file in os.listdir(data_path):
    if curr_file[:-6] in subjects:
        f = open(data_path + curr_file)
        for line in f:
            if line.startswith('>'):
                data.append('>' + line.rstrip().split(' ')[0][1:])
            else:
                if data[-1].startswith('>'):
                    data.append(line.rstrip())
                else:
                    data[-1] += line.rstrip()
        f.close()

f = open(data_path + 'low_data.fasta', 'w')
for line in data:
    if line != '':
        f.write(line + '\n')
f.close()
