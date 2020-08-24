from Source.mtDNA.tibet.functions.file_system import get_path
import pandas as pd

path = get_path()

phylotrees_df = pd.read_excel(path + '/phylotrees.xlsx')
phylotrees_dict = phylotrees_df.to_dict('list')

subjects_info_df = pd.read_excel(path + '/subjects.xlsx')
subjects_info_dict = subjects_info_df.to_dict('list')

fasta_data = {}
classes = ['0-500', '501-1000', '1001-1500', '1501-2000', '2001-2500', '2501-3000', '3001-4000', '4001']
for curr_class in classes:
    f = open(path + '/' + curr_class + '.fasta', 'r')
    curr_subject = ''
    for line in f:
        if line.startswith('>'):
            line_list = line.rstrip().split(' ')
            curr_subject = line_list[0][1:]
            fasta_data[curr_subject] = ''
        else:
            fasta_data[curr_subject] = line.rstrip()
    f.close()

for key in fasta_data:
    curr_dna = fasta_data[key]
    fasta_data[key] = curr_dna.replace('-', '')

positions_to_remove = {}
for key in fasta_data:
    subject_index = subjects_info_dict['subject'].index(key)
    curr_group = subjects_info_dict['group'][subject_index]
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
            curr_index = subjects_info_dict['subject'].index(key)
            if positions_to_remove[position]['group'].startswith(subjects_info_dict['group'][curr_index]):
                curr_dna = fasta_data[key]
                fasta_data[key] = curr_dna[:position] + 'X' + curr_dna[position + 1:]
    elif positions_to_remove[position]['type'] == 'deletion':
        for key in fasta_data:
            curr_index = subjects_info_dict['subject'].index(key)
            if not positions_to_remove[position]['group'].startswith(subjects_info_dict['group'][curr_index]):
                curr_dna = fasta_data[key]
                fasta_data[key] = curr_dna[:position] + 'X' + curr_dna[position + 1:]

fasta_data_mod = {}
for key in fasta_data:
    curr_index = subjects_info_dict['subject'].index(key)
    curr_height = subjects_info_dict['height'][curr_index]
    fasta_data_mod[key + '_' + curr_height] = fasta_data[key].replace('X', '')

f = open(path + '/data_wo_hg.fasta', 'w')
for key in fasta_data_mod:
    f.write('>' + key + '\n')
    f.write(fasta_data_mod[key] + '\n')
f.close()
