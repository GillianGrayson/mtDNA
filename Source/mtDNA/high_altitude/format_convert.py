from Source.mtDNA.tibet.functions.file_system import get_path

path = get_path()

regions_file = 'regions.txt'
regions_info = {}
f = open(path + '/Data/world/rcrs/info/' + regions_file, 'r')
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
fasta_populations = {}
classes = ['Andes', 'Ethiopia', 'Tibetan']
for curr_class in classes:
    fasta_populations[curr_class] = []
    f = open(path + '/Data/world/rcrs/wo_hg/' + curr_class + '.fasta', 'r')
    for line in f:
        if line.startswith('>'):
            line_list = line.rstrip().split(' ')
            fasta_populations[curr_class].append(line_list[0][1:])
            for key in regions_info:
                fasta_data[key].append(line_list[0])
        else:
            for key in regions_info:
                if len(regions_info[key]) > 2:
                    start1 = regions_info[key][0]
                    finish1 = regions_info[key][1] + 1
                    start2 = regions_info[key][2]
                    finish2 = regions_info[key][3] + 1
                    fasta_data[key].append(line.rstrip()[start1:finish1] + line.rstrip()[start2:finish2])
                else:
                    start = regions_info[key][0]
                    finish = regions_info[key][1] + 1
                    fasta_data[key].append(line.rstrip()[start:finish])
    f.close()

for key in fasta_data:
    f = open(path + '/Data/world/rcrs/wo_hg/MT/' + key + '.fasta', 'w')
    for item in fasta_data[key]:
        f.write(item + '\n')
    f.close()

f = open(path + '/Data/world/rcrs/wo_hg/MT/pop.txt', 'w')
for key in fasta_populations:
    for value in fasta_populations[key]:
        f.write(value + '\t' + key + '\n')
f.close()
