import plotly.graph_objs as go
import plotly
import plotly.io as pio
import pandas as pd
import colorlover as cl

def save_figure(fn, fig):
    plotly.offline.plot(fig, filename=fn + '.html', auto_open=False, show_link=True)
    pio.write_image(fig, fn + '.png')
    pio.write_image(fig, fn + '.pdf')

data_path = '../Result/'
data_file_name = 'summary.xlsx'
data = pd.read_excel(data_path + data_file_name)

populations_info = []
for i in range(0, data.shape[0]):
    curr_dict = {}
    curr_dict['long'] = [list(data.longitude)[i]]
    curr_dict['lat'] = [list(data.latitude)[i]]
    curr_dict['name'] = list(data.population_name)[i]
    if list(data.max_precision_of_nuc_therm_snps)[i] == 'reference':
        curr_dict['metrics'] = 'reference'
    else:
        curr_dict['metrics'] = str(list(data.max_precision_of_nuc_therm_snps)[i]) + ' (' + str(list(data.number_of_nuc_therm_snps)[i]) + ')'
    populations_info.append(curr_dict)

color_map = cl.scales['8']['qual']['Set1'] + cl.scales['8']['qual']['Set2']
#color_map = cl.interp(color_map_base, len(populations_info))
colors = color_map[0:len(populations_info)]
coordinates = [color[4:-1].split(',') for color in colors]
colors_transparent = ['rgba(' + ','.join(coordinate) + ',' + str(0.75) + ')' for coordinate in coordinates]

elems = []
for id, population in enumerate(populations_info):
    elems.append(
        go.Scattergeo(
            lon=population['long'],
            lat=population['lat'],
            name=population['name'] + ' - ' + str(population['metrics']),
            hoverinfo='text',
            text=population['name'] + ' - ' + str(population['metrics']),
            mode='markers',
            marker=go.scattergeo.Marker(
                size=5,
                color=colors_transparent[id],
                line=go.scattergeo.marker.Line(
                    width=1,
                    color=colors[id]
                )
            )))

    if id > 0:
        elems.append(
        go.Scattergeo(
            lon = (populations_info[0]['long'][0], population['long'][0]),
            lat = (populations_info[0]['lat'][0], population['lat'][0]),
            name=str(population['metrics']),
            showlegend=False,
            hoverinfo='text',
            text=str(population['metrics']),
            mode = 'lines',
            line = go.scattergeo.Line(
                width = 2,
                color = colors[id],
            )
        ))


layout = go.Layout(
    margin=dict(
        l=10,
        r=10,
        t=10,
        b=10
    ),
    autosize=True,
    title = go.layout.Title(
        text = ''
    ),
    legend=dict(
        x=0.55,
        y=0.49,
        font=dict(
            size=9,
        ),
    ),
    geo = go.layout.Geo(
        showland = True,
        landcolor = 'rgb(243, 243, 243)',
        countrycolor = 'rgb(204, 204, 204)',
        lataxis = go.layout.geo.Lataxis(
            range = [-20, 80],
            dtick = 10
        ),
        lonaxis = go.layout.geo.Lonaxis(
            range = [80, 240],
            dtick = 10
        ),
    ),
)

figure_path = '../Result/figures/'
fig = go.Figure(data=elems , layout=layout)
save_figure(figure_path + 'nuc_therm_KHV_ref', fig)
