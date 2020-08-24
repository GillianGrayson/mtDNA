from Source.mtDNA.tibet.functions.file_system import get_path
from skbio import DistanceMatrix
from skbio.stats.ordination import pcoa
import numpy as np
import pandas as pd

path = get_path()
dm_filename = path + '/FastMe_Distance_Matrix.txt'
dm_array = np.loadtxt(dm_filename)

id_filename = path + '/Mapping.txt'
id_array = np.genfromtxt(id_filename, dtype='str')
ids = [id_array[i][0] for i in range(0, len(id_array))]

metadata = {}
for curr_id in ids:
    curr_list = curr_id.split('_')
    key = curr_list[0]
    if '4001' in curr_list:
        value = '4001'
    else:
        value = curr_list[2] + '-' + curr_list[3]
    metadata[curr_id] = {'height': value}
df = pd.DataFrame.from_dict(metadata, orient='index')

dm = DistanceMatrix(dm_array, ids)
pcoa_results = pcoa(dm)

fig = pcoa_results.plot(df=df, column='height')
fig.show()
