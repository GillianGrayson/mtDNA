import os
import pandas as pd


def get_data(data_path):
    raw_data = []
    subjects = []
    data_classes = []
    for filename in os.listdir(data_path):
        if filename.endswith('fasta'):
            f = open(data_path + filename, 'r')
            raw_data.append([line.rstrip() for line in f][1::2])
            f = open(data_path + filename, 'r')
            subjects.append([line.rstrip().split(' ')[0][1:] for line in f][0::2])
            data_classes.append(filename[:-6])
            f.close()
    return raw_data, subjects, data_classes


def get_haplogroups(data_path):
    haplogroups = {}
    for filename in os.listdir(data_path):
        if filename.endswith('fasta'):
            f = open(data_path + filename, 'r')
            curr_data = [line.rstrip() for line in f][0::2]
            curr_group = filename[:-6]
            haplogroups[curr_group] = []
            for person in curr_data:
                person_list = person.split(' ')
                haplogroup_index = person_list.index('mitochondrion,') - 1
                if 'isolate' in person_list or 'haplogroup' in person_list:
                    haplogroups[curr_group].append(person_list[haplogroup_index])
            f.close()
    return haplogroups


def get_subjects_haplogroups(data_path, subjects):
    haplogroups = {subject: '' for subject_group in subjects for subject in subject_group }
    groups_dict = pd.read_excel(data_path + 'subjects.xlsx').to_dict('list')
    for subject in haplogroups:
        group_index = groups_dict['subject'].index(subject)
        haplogroups[subject] = groups_dict['group'][group_index]
    return haplogroups


def get_regions(data_path):
    regions = {'region': [], 'start': [], 'finish': []}
    f = open(data_path + 'regions.txt', 'r')
    for line in f:
        line_list = line.rstrip().split('\t')
        if line_list[0].startswith('tRNA'):
            line_list[0] = 'tRNA'
        regions['region'].append(line_list[0])
        regions['start'].append(int(line_list[1]))
        regions['finish'].append(int(line_list[2]))
    f.close()
    return regions


def get_mutations_positions(data_path):
    mutation_positions = {}
    phylotree = pd.read_excel(data_path + 'phylotrees.xlsx').to_dict('list')
    for i in range(0, len(phylotree['haplogroup'])):
        if phylotree['haplogroup'][i] not in mutation_positions:
            mutation_positions[phylotree['haplogroup'][i]] = []
        mutation_positions[phylotree['haplogroup'][i]].append(phylotree['position'][i])
    return mutation_positions


def get_features_dict(filename):
    features_dict = {}
    f = open(filename, 'r')
    for line in f:
        line = line.rstrip()
        line_list = line.split('\t')
        features_dict[line_list[0]] = float(line_list[1])
    f.close()
    return features_dict


def save_dict(data, filename):
    f = open(filename, 'w')
    for key in data.keys():
        f.write(str(key) + '\t' + str(data[key]) + '\n')
    f.close()


def save_list(data, filename):
    f = open(filename, 'w')
    for item in data:
        f.write(str(item) + '\n')
    f.close()
