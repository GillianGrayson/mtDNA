from enum import Enum
import socket


class DataPath(Enum):
    local_1 = 'D:/YandexDisk/Work/tibet'
    local_2 = 'E:/YandexDisk/Work/tibet'
    local_3 = 'D:/YandexDisk/tibet'
    local_4 = 'D:/Alena/YandexDisk/tibet'


def get_path():
    host_name = socket.gethostname()
    if host_name == 'MSI':
        path = DataPath.local_1.value
    elif host_name == 'DESKTOP-K9VO2TI':
        path = DataPath.local_2.value
    elif host_name == 'DESKTOP-4BEQ7MS':
        path = DataPath.local_3.value
    else:
        path = DataPath.local_4.value

    return path
