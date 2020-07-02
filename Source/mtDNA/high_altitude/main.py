from Source.mtDNA.tibet.functions.file_system import get_path
from Source.mtDNA.high_altitude.functions import *

path = get_path()
tibet_data_path = path + '/Data/tibet/'
world_data_path = path + '/Data/world/'

tibet_data, tibet_subjects, tibet_classes = read_data(tibet_data_path)
world_data, world_subjects, world_classes = read_data(world_data_path)

regions = get_region_info(tibet_data_path)
olo = 0
