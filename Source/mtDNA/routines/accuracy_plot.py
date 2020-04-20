import os
import json
import hashlib
import plotly.graph_objects as go
import plotly.io as pio
import plotly

data_path = 'E:/YandexDisk/mtDNA/Result/files/'
experiment_type = 'mt-nuc'
random_forest_type = 3
reference_pop = 'IBS'
target_pop = 'TSI'
target_accuracy = 0.5
gene_files = ['mt_gene_list.txt', 'test_gene_list_cold_adaptation.txt']

gene_path = 'E:/YandexDisk/mtDNA/Data/'
genes_mt = []
genes_nuc = []
for file_id in range(0, len(gene_files)):
    data_gene_file = open(gene_path + gene_files[file_id])
    if file_id == 0:
        if experiment_type == 'mt':
            for i, line in enumerate(data_gene_file):
                genes_mt.append(i)
            genes_nuc = []
        elif experiment_type == 'nuc':
            genes_mt = []
            for i, line in enumerate(data_gene_file):
                genes_nuc.append(i)
        elif experiment_type == 'mt-nuc':
            for i, line in enumerate(data_gene_file):
                genes_mt.append(i)
    else:
        for i, line in enumerate(data_gene_file):
            genes_nuc.append(i)
    data_gene_file.close()

json_list = json.dumps([genes_mt, genes_nuc]).encode('utf-8')
curr_hash = hashlib.md5(json_list).hexdigest()

data_path += experiment_type + '/rf_type_' + str(random_forest_type) + \
             '/ref_' + reference_pop + '_target_' + target_pop + '/' + curr_hash + '/'
accuracy_filename = str(target_accuracy) + '_accuracy.txt'
accuracy = []
f = open(data_path + accuracy_filename)
for line in f:
    accuracy.append(float(line.rstrip()))
f.close()

if experiment_type == 'mt' or experiment_type == 'mt-nuc':
    mt_filename = str(target_accuracy) + '_genes_mt.txt'
    mt_genes_num = []
    f = open(data_path + mt_filename)
    for line in f:
        mt_genes_num.append(len(line.rstrip().split('\t')))
    f.close()

if experiment_type == 'nuc' or experiment_type == 'mt-nuc':
    nuc_filename = str(target_accuracy) + '_genes_nuc.txt'
    nuc_genes_num = []
    f = open(data_path + nuc_filename)
    for line in f:
        nuc_genes_num.append(len(line.rstrip().split('\t')))
    f.close()

if experiment_type == 'mt':
    if len(mt_genes_num) > 1:
        if mt_genes_num[0] > mt_genes_num[1]:
            mt_genes_num.append(mt_genes_num.pop(0))
            accuracy.append(accuracy.pop(0))
elif experiment_type == 'nuc':
    if len(nuc_genes_num) > 1:
        if nuc_genes_num[0] > nuc_genes_num[1]:
            nuc_genes_num.append(nuc_genes_num.pop(0))
            accuracy.append(accuracy.pop(0))
else:
    if len(nuc_genes_num) > 1:
        if nuc_genes_num[0] > nuc_genes_num[1]:
            nuc_genes_num.pop(0)
            nuc_genes_num.append(nuc_genes_num[-1] + 1)
            accuracy.append(accuracy.pop(0))

figure_path = 'E:/YandexDisk/mtDNA/Result/figures/'
figure_path += experiment_type + '/rf_type_' + str(random_forest_type) + \
               '/ref_' + reference_pop + '_target_' + target_pop + '/' + curr_hash + '/'

if not os.path.exists(figure_path):
    os.makedirs(figure_path)

fig = go.Figure()
if experiment_type == 'mt':
    fig.add_trace(go.Scatter(x=mt_genes_num, y=accuracy, mode='lines+markers'))
    x_title = 'Number of mtDNA genes'
elif experiment_type == 'nuc':
    fig.add_trace(go.Scatter(x=nuc_genes_num, y=accuracy, mode='lines+markers'))
    x_title = 'Number of nucDNA genes'
else:
    fig.add_trace(go.Scatter(x=nuc_genes_num, y=accuracy, mode='lines+markers'))
    x_title = 'Number of mtDNA-nucDNA gene pairs'

y_title = 'Accuracy'
title = 'Reference ' + reference_pop + ', target ' + target_pop
fig.update_layout(go.Layout(
    title=dict(
        text=title,
        font=dict(
            family='Arial',
            size=25,
        )
    ),
    autosize=True,
    margin=go.layout.Margin(
        l=110,
        r=10,
        b=80,
        t=85,
        pad=0
    ),
    barmode='overlay',
    xaxis=dict(
        title=x_title,
        showgrid=True,
        showline=True,
        mirror='ticks',
        titlefont=dict(
            family='Arial',
            size=25,
            color='black'
        ),
        showticklabels=True,
        tickangle=0,
        tickfont=dict(
            family='Arial',
            size=25,
            color='black'
        )
    ),
    yaxis=dict(
        title=y_title,
        showgrid=True,
        showline=True,
        mirror='ticks',
        titlefont=dict(
            family='Arial',
            size=25,
            color='black'
        ),
        showticklabels=True,
        tickangle=0,
        tickfont=dict(
            family='Arial',
            size=25,
            color='black'
        )
    )
))
plotly.offline.plot(fig, filename=figure_path + 'accuracy.html', auto_open=False, show_link=True)
pio.write_image(fig, figure_path + 'accuracy.png')
pio.write_image(fig, figure_path + 'accuracy.pdf')
