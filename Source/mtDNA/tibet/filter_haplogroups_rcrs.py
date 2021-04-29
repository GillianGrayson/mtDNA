from Source.mtDNA.tibet.functions.file_system import get_path
import pandas as pd

path = get_path()

phylotrees_df = pd.read_excel(path + '/Data/tibet/rcrs/info/phylotrees.xlsx', engine='openpyxl')
phylotrees_dict = phylotrees_df.to_dict('list')

subjects_dict = {'subject': [], 'group': []}

fasta_data = {}
classes = ['0-500', '501-1000', '1001-1500', '1501-2000', '2001-2500', '2501-3000', '3001-4000', '4001']
for curr_class in classes:
    f = open(path + '/Data/tibet/rcrs/rcrs/' + curr_class + '.fasta', 'r')
    curr_subject = ''
    for line in f:
        if line.startswith('>'):
            line_list = line.rstrip().split(' ')
            curr_subject = line_list[0][1:]
            subjects_dict['subject'].append(curr_subject)
            subjects_dict['group'].append(curr_class)
            fasta_data[curr_subject] = ''
        else:
            fasta_data[curr_subject] = line.rstrip()
    f.close()

rcrs_df = pd.read_csv(path + '/Data/tibet/rcrs/info/rCRS_16659.dat', delimiter=',', header=None)
rcrs_dict = rcrs_df.to_dict('list')
rcrs_dict['rcrs'] = rcrs_dict.pop(0)
rcrs_dict['orig'] = rcrs_dict.pop(4)

data = {key: '' for key in fasta_data}
for i in range(0, len(rcrs_dict['orig'])):
    orig_id = int(rcrs_dict['orig'][i])
    rcrs_id = int(rcrs_dict['rcrs'][i])
    if i > 0:
        if rcrs_dict['rcrs'][i] > rcrs_dict['rcrs'][i - 1]:
            for key in fasta_data:
                data[key] += fasta_data[key][orig_id]

f = open(path + '/Data/tibet/rcrs/data_rcrs.fasta', 'w')
for key in data:
    f.write('>' + key + '\n')
    f.write(data[key] + '\n')
f.close()

positions_to_remove = {}
for i in range(0, len(phylotrees_dict['position'])):
    curr_position = phylotrees_dict['position'][i]
    if curr_position not in positions_to_remove:
        if len(str(curr_position).split('-')) == 2:
            for pos in range(int(curr_position.split('-')[0]), int(curr_position.split('-')[1]) + 1):
                positions_to_remove[pos] = {'type': phylotrees_dict['mutation_type'][i],
                                            'ancestral': phylotrees_dict['ancestral_base'][i],
                                            'derived': phylotrees_dict['derived_base'][i],
                                            'group': phylotrees_dict['haplogroup'][i]}
        else:
            positions_to_remove[curr_position] = {'type': phylotrees_dict['mutation_type'][i],
                                                  'ancestral': phylotrees_dict['ancestral_base'][i],
                                                  'derived': phylotrees_dict['derived_base'][i],
                                                  'group': phylotrees_dict['haplogroup'][i]}

positions_list = list(set(positions_to_remove.keys()))
positions_list = [item - 1 for item in positions_list]
positions_list.sort()

data_correspond = {'Original': [], 'New': []}
new_id = 0
for i in range(0, len(data[key])):
    data_correspond['Original'].append(i + 1)
    if (i + 1) not in positions_list:
        new_id += 1
    data_correspond['New'].append(new_id)

info_df = pd.DataFrame(data_correspond)
writer = pd.ExcelWriter(path + '/Data/tibet/rcrs/correspond.xlsx', engine='openpyxl')
info_df.to_excel(writer, index=False, startrow=0)
worksheet = writer.sheets['Sheet1']
writer.save()

data_wo_hg = {}
for key in data:
    curr_dna = data[key]
    data_wo_hg[key] = ''.join([i for j, i in enumerate(curr_dna) if j not in positions_list])

f = open(path + '/Data/tibet/rcrs/data_wo_hg.fasta', 'w')
for key in data_wo_hg:
    f.write('>' + key + '\n')
    f.write(data_wo_hg[key] + '\n')
f.close()

for curr_class in classes:
    f = open(path + '/Data/tibet/rcrs/wo_hg/' + curr_class + '.fasta', 'w')
    curr_ids = [i for i, x in enumerate(subjects_dict['group']) if x == curr_class]
    for curr_id in curr_ids:
        f.write('>' + subjects_dict['subject'][curr_id] + '\n')
        f.write(data_wo_hg[subjects_dict['subject'][curr_id]] + '\n')
    f.close()
