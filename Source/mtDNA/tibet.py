import os
import itertools
from Source.mtDNA.tibet_functions import create_df_freq, create_df, run_random_forest

use_freq = 1

data_path = 'C:/Users/User/YandexDisk/tibet/Data/'
result_path = 'C:/Users/User/YandexDisk/tibet/Result/'

raw_data = []
data_classes = []
for filename in os.listdir(data_path):
    f = open(data_path + filename, 'r')
    raw_data.append([line.rstrip() for line in f][1::2])
    data_classes.append(filename[:-6])
    f.close()

if use_freq:
    df_main = create_df_freq(raw_data)
else:
    df_main = create_df(raw_data)

classes = []
for i in range(0, len(raw_data)):
    classes += [data_classes[i], ] * len(raw_data[i])

accuracy = run_random_forest(df_main, classes)
print('8-class Classification Accuracy: ' + str(accuracy))

result_file_name = 'classification.txt'
if not os.path.exists(result_path):
    os.makedirs(result_path)
f = open(result_path + result_file_name, 'w')
f.write('8-class Classification Accuracy: ' + str(accuracy) + '\n')

for subset in itertools.combinations(raw_data, 2):
    test_data = subset
    if use_freq:
        df_main = create_df_freq(test_data)
    else:
        df_main = create_df(test_data)

    classes = [data_classes[raw_data.index(test_data[0])], ] * len(test_data[0]) + \
              [data_classes[raw_data.index(test_data[1])], ] * len(test_data[1])

    accuracy = run_random_forest(df_main, classes)
    print(data_classes[raw_data.index(test_data[0])] + ' vs ' + data_classes[raw_data.index(test_data[1])] +
          ' Binary Classification Accuracy: ' + str(accuracy))
    f.write(data_classes[raw_data.index(test_data[0])] + ' vs ' + data_classes[raw_data.index(test_data[1])] +
            ' Binary Classification Accuracy: ' + str(accuracy) + '\n')
f.close()
