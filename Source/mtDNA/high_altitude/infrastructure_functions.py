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


def save_results(path, filename, data):
    f = open(path + filename + '.txt', 'w')
    for item in data:
        f.write(str(item) + '\n')
    f.close()


def read_results(path, filename):
    data = []
    f = open(path + filename, 'r')
    for line in f:
        line = line.rstrip()
        data.append(float(line))
    f.close()
    return data
