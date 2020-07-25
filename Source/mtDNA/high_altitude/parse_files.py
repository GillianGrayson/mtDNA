from Source.mtDNA.tibet.functions.file_system import get_path

path = get_path()
aligned_data_path = path + '/Data/alignment/'

tibet_classes = ['0-500', '501-1000', '1001-1500', '1501-2000', '2001-2500', '2501-3000', '3001-4000', '4001']
world_classes = ['Andes', 'Ethiopia', 'Tibetan']

data = []
for tibet_class in tibet_classes:
    f = open(aligned_data_path + tibet_class + '.fasta')
    for line in f:
        data.append(line.rstrip())
    f.close()

for world_class in world_classes:
    f = open(aligned_data_path + world_class + '.fa')
    subject_count = 1
    for line in f:
        if line.startswith('>'):
            data.append('>' + world_class + ' ' + str(subject_count))
            subject_count += 1
        else:
            data.append(line.rstrip())
    f.close()

f = open(aligned_data_path + 'all_data.fasta', 'w')
for line in data:
    if line != '':
        f.write(line + '\n')
f.close()
