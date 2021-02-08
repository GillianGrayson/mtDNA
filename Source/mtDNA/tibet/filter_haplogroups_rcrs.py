from Source.mtDNA.tibet.functions.file_system import get_path
import pandas as pd

path = get_path()

phylotrees_df = pd.read_excel(path + '/info/phylotrees.xlsx', engine='openpyxl')
phylotrees_dict = phylotrees_df.to_dict('list')

subjects_dict = {'subject': [], 'group': []}

fasta_data = {}
classes = ['Andes', 'Ethiopia', 'Tibetan']
for curr_class in classes:
    f = open(path + '/' + curr_class + '.fasta', 'r')
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

rcrs_df = pd.read_csv(path + '/info/rCRS_TEA.dat', delimiter='\t', header=None)
rcrs_dict = rcrs_df.to_dict('list')
rcrs_dict['rcrs'] = rcrs_dict.pop(0)
rcrs_dict['orig'] = rcrs_dict.pop(1)

data = {key: '' for key in fasta_data}
for i in range(0, len(rcrs_dict['orig'])):
    orig_id = rcrs_dict['orig'][i]
    rcrs_id = rcrs_dict['rcrs'][i]
    if i > 0:
        if rcrs_dict['rcrs'][i] > rcrs_dict['rcrs'][i - 1]:
            for key in fasta_data:
                data[key] += fasta_data[key][orig_id]

f = open(path + '/data_rcrs.fasta', 'w')
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
data_wo_hg = {}
for key in data:
    curr_dna = data[key]
    data_wo_hg[key] = ''.join([i for j, i in enumerate(curr_dna) if j not in positions_list])

f = open(path + '/data_wo_hg.fasta', 'w')
for key in data_wo_hg:
    f.write('>' + key + '\n')
    f.write(data_wo_hg[key] + '\n')
f.close()
