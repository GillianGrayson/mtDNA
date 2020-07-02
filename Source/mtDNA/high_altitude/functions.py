import os


def read_data(data_path):
    raw_data = []
    subjects = []
    data_classes = []
    for filename in os.listdir(data_path):
        if filename.endswith('fasta') or filename.endswith('fa'):
            f = open(data_path + filename, 'r')
            raw_data.append([line.rstrip() for line in f][1::2])
            f = open(data_path + filename, 'r')
            subjects.append([line.rstrip().split(' ')[0][1:] for line in f][0::2])
            if filename.endswith('fasta'):
                data_classes.append(filename[:-6])
            else:
                data_classes.append(filename[:-3])
            f.close()
    return raw_data, subjects, data_classes


def get_region_info(data_path):
    regions = {'region': [], 'start': [], 'finish': []}
    f = open(data_path + 'regions.txt', 'r')
    for line in f:
        line_list = line.rstrip().split('\t')
        regions['region'].append(line_list[0])
        regions['start'].append(int(line_list[1]))
        regions['finish'].append(int(line_list[2]))
    f.close()
    return regions
