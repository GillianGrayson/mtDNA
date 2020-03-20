from Source.mtDNA.tibet.functions.file_system import get_path
import os

path = get_path()
data_path = path + '/Data/'
result_path = path + '/Data/haplogroups/'

if not os.path.exists(result_path):
    os.makedirs(result_path)

haplogroups_classes = []
f = open(data_path + 'haplogroups.txt')
for line in f:
    haplogroups_classes.append([item for item in line.rstrip().split(',')])
f.close()

raw_data = []
raw_haplogroups = []
raw_persons = []
for filename in os.listdir(data_path):
    if filename.endswith('fasta'):
        f = open(data_path + filename, 'r')
        for line in f:
            if line.startswith('>'):
                person_list = line.split(' ')
                haplogroup_index = person_list.index('mitochondrion,') - 1
                if 'isolate' in person_list or 'haplogroup' in person_list:
                    raw_haplogroups.append(person_list[haplogroup_index])
                    raw_persons.append(person_list[0])
                    is_haplogroup = 1
                else:
                    is_haplogroup = 0
            else:
                if is_haplogroup:
                    raw_data.append(line.rstrip())
        f.close()

for haplogroups_class in haplogroups_classes:
    if '*' in haplogroups_class[0]:
        suffix = haplogroups_class[0].replace('*', '')
    else:
        suffix = haplogroups_class[0]
    f = open(result_path + suffix + '.fasta', 'w')
    for haplogroup in haplogroups_class:
        index = raw_haplogroups.index(haplogroup)
        f.write(raw_persons[index] + ' ' + raw_haplogroups[index] + '\n')
        f.write(raw_data[index] + '\n')
    f.close()
