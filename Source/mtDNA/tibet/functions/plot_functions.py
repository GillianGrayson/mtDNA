import plotly
import plotly.graph_objs as go
import pickle


def plot_hist(data, suffix, file_path):
    fig = go.Figure(go.Bar(
        x=list(data.values())[::-1],
        y=list(data.keys())[::-1],
        orientation='h'
    ))
    fig.update_yaxes(
        tickfont=dict(size=10)
    )
    fig.update_layout(width=700,
                      height=1000)

    plotly.offline.plot(fig, filename=file_path + suffix + '_hist.html', auto_open=False, show_link=True)
    # plotly.io.write_image(fig, file_path + suffix + '_hist.png')
    # plotly.io.write_image(fig, file_path + suffix + '_hist.pdf')

    with open(file_path + suffix + '_hist_fig.pickle', 'wb') as handle:
        pickle.dump(fig, handle, protocol=pickle.HIGHEST_PROTOCOL)


def plot_hist_from_pkl(suffix, file_path):
    with open(file_path + suffix + '_hist_fig.pickle', 'rb') as handle:
        fig_pkl = pickle.load(handle)

    plotly.offline.plot(fig_pkl, filename=file_path + suffix + '_hist.html', auto_open=False, show_link=True)
    plotly.io.write_image(fig_pkl, file_path + suffix + '_hist.png')
    plotly.io.write_image(fig_pkl, file_path + suffix + '_hist.pdf')
