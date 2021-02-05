from Source.mtDNA.tibet.functions.file_system import get_path
import pandas as pd
from collections import Counter

path = get_path()

phylotrees_df = pd.read_excel(path + '/phylotrees.xlsx')
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

rcrs_df = pd.read_csv(path + '/rCRS_wo_gap.dat', delimiter='\t', header=None)
rcrs_dict = rcrs_df.to_dict('list')
rcrs_dict['rcrs'] = rcrs_dict.pop(0)
rcrs_dict['orig'] = rcrs_dict.pop(1)
rcrs_dict['nuc'] = rcrs_dict.pop(2)

curr_rcrs_id = 0
ids_to_remove = []
for i in range(0, len(fasta_data[subjects_dict['subject'][0]])):
    curr_nucs = []
    for key in fasta_data:
        curr_nucs.append(fasta_data[key][i])
    count_dict = Counter(curr_nucs).most_common()
    count_dict = dict(count_dict)
    nucs = list(count_dict.keys())
    if rcrs_dict['nuc'][curr_rcrs_id] in nucs and curr_rcrs_id + 1 == rcrs_dict['rcrs'][curr_rcrs_id]:
        curr_rcrs_id += 1
    else:
        ids_to_remove.append(i)

for key in fasta_data:
    curr_dna = fasta_data[key]
    fasta_data[key] = ''.join([i for j, i in enumerate(curr_dna) if j not in ids_to_remove])

for key in fasta_data:
    curr_dna = fasta_data[key]
    fasta_data[key] = curr_dna.replace('-', '')

positions_to_remove = {}
for key in fasta_data:
    subject_index = subjects_dict['subject'].index(key)
    curr_group = subjects_dict['group'][subject_index]
    precise = True
    if curr_group.endswith('*'):
        curr_group = curr_group[:-1]
        precise = False
    group_indices = []
    for i in range(0, len(phylotrees_dict['haplogroup'])):
        group = phylotrees_dict['haplogroup'][i]
        if precise:
            if group == curr_group:
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
        else:
            if group.startswith(curr_group):
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

positions_list = list(positions_to_remove.keys())
positions_list.sort()
for position in positions_list:
    if positions_to_remove[position]['type'] == 'mutation':
        for key in fasta_data:
            curr_dna = fasta_data[key]
            fasta_data[key] = curr_dna[:position] + 'X' + curr_dna[position+1:]
    elif positions_to_remove[position]['type'] == 'insertion':
        for key in fasta_data:
            curr_index = subjects_dict['subject'].index(key)
            if positions_to_remove[position]['group'].startswith(subjects_dict['group'][curr_index]):
                curr_dna = fasta_data[key]
                fasta_data[key] = curr_dna[:position] + 'X' + curr_dna[position + 1:]
    elif positions_to_remove[position]['type'] == 'deletion':
        for key in fasta_data:
            curr_index = subjects_dict['subject'].index(key)
            if not positions_to_remove[position]['group'].startswith(subjects_dict['group'][curr_index]):
                curr_dna = fasta_data[key]
                fasta_data[key] = curr_dna[:position] + 'X' + curr_dna[position + 1:]

fasta_data_mod = {}
for key in fasta_data:
    curr_index = subjects_dict['subject'].index(key)
    fasta_data_mod[key] = fasta_data[key].replace('X', '')

f = open(path + '/data_wo_hg.fasta', 'w')
for key in fasta_data_mod:
    f.write('>' + key + '\n')
    f.write(fasta_data_mod[key] + '\n')
f.close()
