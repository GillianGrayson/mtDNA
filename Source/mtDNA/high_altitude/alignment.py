from Source.mtDNA.tibet.functions.file_system import get_path
import os

read_tibet_data = 1

path = get_path()
tibet_data_path = path + '/Data/tibet/'
world_data_path = path + '/Data/world/'
aligned_data_path = path + '/Data/alignment/'

tibet_subjects_count = {'0-500': 0, '501-1000': 0, '1001-1500': 0, '1501-2000': 0, '2001-2500': 0, '2501-3000': 0,
                        '3001-4000': 0, '4001': 0}

for filename in os.listdir(tibet_data_path):
    if filename.endswith('fasta'):
        f = open(tibet_data_path + filename, 'r')
        for i, l in enumerate(f):
            pass
        f.close()
        tibet_subjects_count[filename[:-6]] = (i + 1) // 2

world_subjects_count = {'Andes': 0, 'Ethiopia': 0, 'Tibetan': 0}

for filename in os.listdir(world_data_path):
    if filename.endswith('fa'):
        f = open(world_data_path + filename, 'r')
        for i, l in enumerate(f):
            pass
        f.close()
        world_subjects_count[filename[:-3]] = (i + 1) // 2

data = []
f = open(aligned_data_path + 'all_data_wo_gaps_aligned.fasta')
for line in f:
    if line.startswith('>'):
        data.append(line.rstrip())
        prev_line = line.rstrip()
    else:
        if prev_line.startswith('>'):
            data.append(line.rstrip())
            prev_line = line.rstrip()
        else:
            data[-1] += line.rstrip()
            prev_line = line.rstrip()
f.close()

aligned_data = {'tibet': {}, 'world': {}}
tibet_groups = list(tibet_subjects_count.keys())
global_count = 0
num_lines = 0
for tibet_group in tibet_groups:
    aligned_data['tibet'][tibet_group] = []
    num_lines += tibet_subjects_count[tibet_group] * 2
    while global_count < num_lines:
        aligned_data['tibet'][tibet_group].append(data[global_count])
        global_count += 1

world_groups = list(world_subjects_count.keys())
for world_group in world_groups:
    aligned_data['world'][world_group] = []
    num_lines += world_subjects_count[world_group] * 2
    while global_count < num_lines:
        aligned_data['world'][world_group].append(data[global_count])
        global_count += 1

tibet_result_path = aligned_data_path + 'tibet/'
if not os.path.exists(tibet_result_path):
    os.makedirs(tibet_result_path)
for group in aligned_data['tibet']:
    f = open(tibet_result_path + group + '.fasta', 'w')
    for item in aligned_data['tibet'][group]:
        f.write(str(item) + '\n')
    f.close()

world_result_path = aligned_data_path + 'world/'
if not os.path.exists(world_result_path):
    os.makedirs(world_result_path)
for group in aligned_data['world']:
    f = open(world_result_path + group + '.fasta', 'w')
    for item in aligned_data['world'][group]:
        f.write(str(item) + '\n')
    f.close()
