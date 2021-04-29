import os
import xlsxwriter
import pandas as pd


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
        data.append(line)
    f.close()
    return data


def write_frequency_to_xlsx(path, filename, data):
    workbook = xlsxwriter.Workbook(path + filename + '.xlsx')
    worksheet = workbook.add_worksheet()
    row_id = 0
    column_id = 0
    if len(data) > 0:
        keys = list(data.keys())
        column_names = list(data[keys[0]].keys())
        worksheet.write(row_id, column_id, 'Id')
        for key_id in range(0, len(column_names)):
            worksheet.write(row_id, column_id + key_id + 1, column_names[key_id])
        row_id += 1
        for key in data:
            worksheet.write(row_id, column_id, key)
            for column_name_id in range(0, len(column_names)):
                try:
                    worksheet.write(row_id, column_id + column_name_id + 1, data[key][column_names[column_name_id]])
                except:
                    pass
            row_id += 1
    workbook.close()


def write_regions_to_xlsx(path, filename, data):
    workbook = xlsxwriter.Workbook(path + filename + '.xlsx')
    worksheet = workbook.add_worksheet()
    row_id = 0
    column_id = 0
    worksheet.write(row_id, column_id, 'Region')
    worksheet.write(row_id, column_id + 1, 'Frequency')
    row_id += 1
    if len(data) > 0:
        for key in data:
            worksheet.write(row_id, column_id, key)
            worksheet.write(row_id, column_id + 1, data[key])
            row_id += 1
    workbook.close()


def write_stat_to_xlsx(path, filename, data, num_lines):
    df = pd.DataFrame(data)
    if num_lines == 'all':
        new_df = df
    elif isinstance(num_lines, int):
        new_df = df.head(num_lines)
    new_df.to_excel(path + filename + '.xlsx', index=False)
