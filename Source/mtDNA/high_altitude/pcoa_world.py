from Source.mtDNA.tibet.functions.file_system import get_path
from skbio import DistanceMatrix
from skbio.stats.ordination import pcoa
import numpy as np
import pandas as pd
import plotly
# import colorlover as cl
import plotly.graph_objs as go
# import plotly.express as px
import cmocean


def cmocean_to_plotly(cmap, pl_entries):
    h = 1.0/(pl_entries-1)
    pl_colorscale = []

    for k in range(pl_entries):
        C = list(map(np.uint8, np.array(cmap(k*h)[:3])*255))
        pl_colorscale.append([k*h, 'rgb'+str((C[0], C[1], C[2]))])

    return pl_colorscale


path = get_path()
dm_filename = path + '/DM.txt'
dm_array = np.loadtxt(dm_filename)

id_filename = path + '/Mapping.txt'
id_array = np.genfromtxt(id_filename, dtype='str')
ids = [id_array[i][0] for i in range(0, len(id_array))]

metadata = {}
subject_dict = {}
for curr_id in ids:
    curr_list = curr_id.split('_')
    key = curr_list[0]
    if key in subject_dict:
        subject_dict[key].append(curr_id)
    else:
        subject_dict[key] = [curr_id]
    metadata[curr_id] = {'group': key}
df = pd.DataFrame.from_dict(metadata, orient='index')

dm = DistanceMatrix(dm_array, ids)
pcoa_results = pcoa(dm)

fig = pcoa_results.plot(df=df, column='group')
fig.show()

coord_matrix = pcoa_results.samples.values.T
xs_all = coord_matrix[0]
ys_all = coord_matrix[1]

traces_2d = []
for status in subject_dict:
    curr_subjects = subject_dict[status]
    xs = []
    ys = []
    for subj in curr_subjects:
        index = ids.index(subj)
        xs.append(xs_all[index])
        ys.append(ys_all[index])

    color = cmocean_to_plotly(cmocean.cm.phase, 4)[list(subject_dict.keys()).index(status)]
    coordinates = color[1][4:-1].split(',')
    color_transparent = 'rgba(' + ','.join(['0', '0', '0']) + ',' + str(0.7) + ')'
    color_border = 'rgba(' + ','.join(coordinates) + ',' + str(1) + ')'

    # color = cl.scales['8']['qual']['Set1'][list(subject_dict.keys()).index(status)]
    # coordinates = color[4:-1].split(',')
    # color_transparent = 'rgba(' + ','.join(coordinates) + ',' + str(0.3) + ')'
    # color_border = 'rgba(' + ','.join(coordinates) + ',' + str(1) + ')'

    trace = go.Scatter(
        x=ys,
        y=xs,
        name=status,
        mode='markers',
        marker=dict(
            size=12,
            color=color_border,
            line=dict(
                color=color_transparent,
                width=1.0
            ),
            opacity=1.0
        )
    )
    traces_2d.append(trace)

layout_2d = go.Layout(
    plot_bgcolor='rgba(0,0,0,0)',
    margin=go.layout.Margin(
        l=0,
        r=0,
        b=0,
        t=0,
        pad=0
    ),
    autosize=True,
    legend=dict(
        font=dict(
            family='Arial',
            size=16,
        )
    ),
    xaxis=dict(
        title='PC1',
        showgrid=True,
        showline=True,
        linewidth=1,
        linecolor='rgba(0,0,0,0.5)',
        gridcolor='rgba(0,0,0,0.1)',
        mirror=True,
        titlefont=dict(
            family='Arial',
            color='black',
            size=24,
        ),
        showticklabels=True,
        tickangle=0,
        tickfont=dict(
            family='Arial',
            color='black',
            size=20
        ),
        exponentformat='e',
        showexponent='all',
        zeroline=True,
        zerolinecolor='rgba(0,0,0,0.1)'
    ),
    yaxis=dict(
        title='PC2',
        showgrid=True,
        showline=True,
        linewidth=1,
        linecolor='rgba(0,0,0,0.5)',
        gridcolor='rgba(0,0,0,0.1)',
        mirror=True,
        titlefont=dict(
            family='Arial',
            color='black',
            size=24,
        ),
        showticklabels=True,
        tickangle=0,
        tickfont=dict(
            family='Arial',
            color='black',
            size=20
        ),
        exponentformat='e',
        showexponent='all',
        zeroline=True,
        zerolinecolor='rgba(0,0,0,0.1)'
    ),
)

fig_2d = go.Figure(data=traces_2d, layout=layout_2d)

plotly.offline.plot(fig_2d, filename=path + '/pcoa_2d.html', auto_open=False, show_link=True)
plotly.io.write_image(fig_2d, path + '/pcoa_2d.png')
plotly.io.write_image(fig_2d, path + '/pcoa_2d.pdf')
