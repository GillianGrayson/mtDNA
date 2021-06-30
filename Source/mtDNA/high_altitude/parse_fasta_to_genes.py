from Source.mtDNA.tibet.functions.file_system import get_path
import pandas as pd


path = get_path()
alignment_data_path = path + '/Data/alignment/low_data/'
info_data_path = path + '/Data/alignment/info/'
result_data_path = path + '/Result/low_data/genes/'

data_dict = {}
f = open(alignment_data_path + 'all_low_wo_rcrs.txt', 'r')
for line in f:
    if line.startswith('>'):
        curr_key = line.rstrip().split(' ')[0][1:]
        data_dict[curr_key] = ''
    elif prev_line.startswith('>'):
        data_dict[curr_key] = line.rstrip()
    else:
        data_dict[curr_key] += line.rstrip()
    prev_line = line
f.close()

subjects_data = pd.read_excel(info_data_path + 'subjects_all.xlsx').to_dict('list')
high_classes = ['Andes', 'Tibetan', 'Ethiopia']
subject_info = {high_class: [] for high_class in high_classes}
for subject in data_dict:
    if subject.split('_')[0] in high_classes:
        subject_info[subject.split('_')[0]].append(subject)
    else:
        subject_id = subjects_data['subject'].index(subject)
        subject_group = subjects_data['height'][subject_id]
        if subject_group not in subject_info:
            subject_info[subject_group] = [subject]
        else:
            subject_info[subject_group].append(subject)

regions_file = 'regions.txt'
regions_info = {}
f = open(info_data_path + regions_file, 'r')
for line in f:
    line_list = line.rstrip().split('\t')
    if line_list[0] in regions_info:
        regions_info[line_list[0]].append(int(line_list[1]))
        regions_info[line_list[0]].append(int(line_list[2]))
    else:
        regions_info[line_list[0]] = [int(line_list[1]), int(line_list[2])]
f.close()

fasta_data = {}
fasta_data.update({key: [] for key in list(regions_info.keys())})
for subject in data_dict:
    for region in regions_info:
        fasta_data[region].append('>' + subject)
        if len(regions_info[region]) > 2:
            start1 = regions_info[region][0]
            finish1 = regions_info[region][1] + 1
            start2 = regions_info[region][2]
            finish2 = regions_info[region][3] + 1
            fasta_data[region].append(data_dict[subject][start1:finish1] + data_dict[subject][start2:finish2])
        else:
            start = regions_info[region][0]
            finish = regions_info[region][1] + 1
            fasta_data[region].append(data_dict[subject][start:finish])

for key in fasta_data:
    f = open(result_data_path + key + '.fasta', 'w')
    for item in fasta_data[key]:
        f.write(item + '\n')
    f.close()

f = open(result_data_path + 'pop.txt', 'w')
for group in subject_info:
    for subject in subject_info[group]:
        if group == '0-500':
            group_val = 'Low'
        else:
            group_val = group
        f.write(subject + '\t' + group_val + '\n')
f.close()
