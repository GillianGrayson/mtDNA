import plotly
import plotly.io as pio
from plotly.graph_objs import *


def save_figure(fn, fig):
    plotly.offline.plot(fig, filename=fn + '.html', auto_open=False, show_link=True)
    pio.write_image(fig, fn + '.png')
    pio.write_image(fig, fn + '.pdf')


result_path = 'E:/YandexDisk/mtDNA/Result/figures/map/'

trace1 = {
    "type": "choropleth",
    "showscale": False,
    "locationmode": "country names",
    "locations": ['Great Britain', 'Finland', 'Italy'],
    "z": [0.0, 0.5, 1.0],
    "text": ['GBR', 'FIN', 'TSI'],
    "colorscale": 'Viridis'
}
layout = {
    "geo": {
        "scope": "europe",
        "domain": {
            "x": [0, 1],
            "y": [0, 1]
        },
        "lataxis": {"range": [35.0, 75.0]},
        "lonaxis": {"range": [-15.0, 35.0]},
        "showland": True,
        "landcolor": "rgb(229, 229, 229)",
        "showframe": True,
        "projection": {"type": "mercator"},
        "resolution": 50,
        "countrycolor": "rgb(255, 0, 255)",
        "coastlinecolor": "rgb(0, 255, 255)",
        "showcoastlines": True,
        "showcountries": True
    },
    "legend": {"traceorder": "reversed"}
}
fig = Figure(data=[trace1], layout=layout)
save_figure(result_path + 'Europe', fig)
