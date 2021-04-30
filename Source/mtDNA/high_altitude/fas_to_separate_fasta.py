from Source.mtDNA.tibet.functions.file_system import get_path
import os

path = get_path()
aligned_data_path = path + '/Data/alignment/align_deletion/'
aligned_file = path + '/Data/alignment/align_deletion.fas'

data = []
f = open(aligned_file, 'r')
for line in f:
    if line.startswith('>'):
        data.append(line.rstrip())
    elif data[-1].startswith('>'):
        data.append(line.rstrip())
    else:
        data[-1] += line.rstrip()
f.close()

tibet_subjects_count = {'0-500': 0, '501-1000': 0, '1001-1500': 0, '1501-2000': 0, '2001-2500': 0, '2501-3000': 0,
                        '3001-4000': 0, '4001': 0}

for filename in os.listdir(path + '/Data/alignment/tibet/'):
    if filename.endswith('fasta'):
        f = open(path + '/Data/alignment/tibet/' + filename, 'r')
        for i, l in enumerate(f):
            pass
        f.close()
        tibet_subjects_count[filename[:-6]] = (i + 1) // 2

world_subjects_count = {'Andes': 0, 'Ethiopia': 0, 'Tibetan': 0}

for filename in os.listdir(path + '/Data/alignment/world/'):
    if filename.endswith('fasta'):
        f = open(path + '/Data/alignment/world/' + filename, 'r')
        for i, l in enumerate(f):
            pass
        f.close()
        world_subjects_count[filename[:-6]] = (i + 1) // 2

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
    f = open(tibet_result_path + group + '.fas', 'w')
    for item in aligned_data['tibet'][group]:
        f.write(str(item) + '\n')
    f.close()

world_result_path = aligned_data_path + 'world/'
if not os.path.exists(world_result_path):
    os.makedirs(world_result_path)
for group in aligned_data['world']:
    f = open(world_result_path + group + '.fas', 'w')
    for item in aligned_data['world'][group]:
        f.write(str(item) + '\n')
    f.close()
