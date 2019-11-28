import os


def get_data(data_path):
    raw_data = []
    data_classes = []
    for filename in os.listdir(data_path):
        if filename.endswith('fasta'):
            f = open(data_path + filename, 'r')
            raw_data.append([line.rstrip() for line in f][1::2])
            data_classes.append(filename[:-6])
            f.close()
    return raw_data, data_classes


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
