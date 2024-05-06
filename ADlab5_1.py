import dash
from dash import dcc, html, Input, Output
import numpy as np

app = dash.Dash(__name__)

initial_amplitude = 1.0
initial_frequency = 1.0
initial_phase = 0.0
initial_noise_mean = 0.0
initial_noise_covariance = 0.1
initial_window_size = 5

# Keeping track of previous noise and its parameters
previous_noise = None
previous_noise_mean = None
previous_noise_covariance = None

t = np.arange(0.0, 10.0, 0.01)

def harmonic_with_noise(t, amplitude, frequency, phase, noise_mean, noise_covariance, show_noise):
    signal = amplitude * np.sin(2 * np.pi * frequency * t + phase)
    if show_noise:
        noise = np.random.normal(noise_mean, np.sqrt(noise_covariance), len(t))
        return signal + noise, signal
    else:
        return signal, signal

def custom_filter(signal, window_size):
    filtered_signal = np.zeros_like(signal)
    half_window = window_size // 2

    for i in range(len(signal)):
        start = max(0, i - half_window)
        end = min(len(signal), i + half_window + 1)
        window_sum = 0
        window_length = end - start

        for j in range(start, end):
            window_sum += signal[j]

        filtered_signal[i] = window_sum / window_length

    return filtered_signal

y_noise, y_harmonic = harmonic_with_noise(t, initial_amplitude, initial_frequency, initial_phase, initial_noise_mean, initial_noise_covariance, True)
filtered_signal = custom_filter(y_noise, initial_window_size)

app.layout = html.Div([
    html.H1("Harmonic Signal with Noise"),
    dcc.Graph(id='main-graph', figure={
        'data': [
            {'x': t, 'y': y_noise, 'type': 'scatter', 'name': 'Signal with Noise'},
            {'x': t, 'y': y_harmonic, 'type': 'scatter', 'name': 'Pure Harmonic'},
            {'x': t, 'y': filtered_signal, 'type': 'scatter', 'name': 'Filtered Signal'}
        ],
        'layout': {
            'title': 'Harmonic Signal with Noise',
            'xaxis': {'title': 'Time'},
            'yaxis': {'title': 'Amplitude'}
        }
    }),
    html.Label("Amplitude"),
    dcc.Slider(id='amplitude-slider', min=0.1, max=2.0, step=0.1, value=initial_amplitude),
    html.Label("Frequency"),
    dcc.Slider(id='frequency-slider', min=0.1, max=2.0, step=0.1, value=initial_frequency),
    html.Label("Phase"),
    dcc.Slider(id='phase-slider', min=0, max=2*np.pi, step=0.1, value=initial_phase),
    html.Label("Noise Mean"),
    dcc.Slider(id='noise-mean-slider', min=-1.0, max=1.0, step=0.1, value=initial_noise_mean),
    html.Label("Noise Covariance"),
    dcc.Slider(id='noise-covariance-slider', min=0.0, max=1.0, step=0.1, value=initial_noise_covariance),
    html.Label("Window Size"),
    dcc.Slider(id='window-size-slider', min=0, max=100, step=5, value=initial_window_size),
    html.Button('Reset', id='reset-button', n_clicks=0),
    html.Button('Show Noise', id='show-noise-button', n_clicks=0),
    html.Label("Select Graph to Display"),
    dcc.Dropdown(
        id='graph-selector',
        options=[
            {'label': 'Signal with Noise', 'value': 'signal_noise'},
            {'label': 'Pure Harmonic', 'value': 'pure_harmonic'},
            {'label': 'Filtered Signal', 'value': 'filtered_signal'}
        ],
        value='signal_noise'
    ),
    html.Div([
        dcc.Graph(id='secondary-graph')
    ])
])

@app.callback(
    [
        Output('amplitude-slider', 'value'),
        Output('frequency-slider', 'value'),
        Output('phase-slider', 'value'),
        Output('noise-mean-slider', 'value'),
        Output('noise-covariance-slider', 'value'),
        Output('window-size-slider', 'value')
    ],
    [Input('reset-button', 'n_clicks')]
)
def reset_sliders(reset_clicks):
    if reset_clicks:
        return initial_amplitude, initial_frequency, initial_phase, initial_noise_mean, initial_noise_covariance, initial_window_size
    else:
        raise dash.exceptions.PreventUpdate

@app.callback(
    Output('main-graph', 'figure'),
    [
        Input('amplitude-slider', 'value'),
        Input('frequency-slider', 'value'),
        Input('phase-slider', 'value'),
        Input('noise-mean-slider', 'value'),
        Input('noise-covariance-slider', 'value'),
        Input('window-size-slider', 'value'),
        Input('show-noise-button', 'n_clicks')
    ],
    prevent_initial_call=True
)
def update_graph(amplitude, frequency, phase, noise_mean, noise_covariance, window_size, noise_clicks):
    global previous_noise, previous_noise_mean, previous_noise_covariance

    if (noise_mean != previous_noise_mean or noise_covariance != previous_noise_covariance):
        previous_noise_mean = noise_mean
        previous_noise_covariance = noise_covariance
        previous_noise = np.random.normal(noise_mean, np.sqrt(noise_covariance), len(t))

    y_noise, y_harmonic = amplitude * np.sin(2 * np.pi * frequency * t + phase), amplitude * np.sin(2 * np.pi * frequency * t + phase)
    
    show_noise = noise_clicks % 2 == 1

    if show_noise:
        y_noise += previous_noise
        filtered_signal = custom_filter(y_noise, window_size)
    else:
        filtered_signal = None

    fig = {
        'data': [
            {'x': t, 'y': y_noise, 'type': 'scatter', 'name': 'Signal with Noise'},
            {'x': t, 'y': y_harmonic, 'type': 'scatter', 'name': 'Pure Harmonic'}
        ],
        'layout': {
            'title': 'Harmonic Signal with Noise',
            'xaxis': {'title': 'Time'},
            'yaxis': {'title': 'Amplitude'}
        }
    }

    if filtered_signal is not None:
        fig['data'].append({'x': t, 'y': filtered_signal, 'type': 'scatter', 'name': 'Filtered Signal'})

    return fig

@app.callback(
    Output('secondary-graph', 'figure'),
    [
        Input('graph-selector', 'value'),
        Input('main-graph', 'figure')
    ],
    prevent_initial_call=True
)
def update_secondary_graph(selected_graph, main_graph_figure):
    main_data = main_graph_figure['data']
    if selected_graph == 'signal_noise':
        return {'data': [trace for trace in main_data if trace['name'] == 'Signal with Noise']}
    elif selected_graph == 'pure_harmonic':
        return {'data': [trace for trace in main_data if trace['name'] == 'Pure Harmonic']}
    elif selected_graph == 'filtered_signal':
        return {'data': [trace for trace in main_data if trace['name'] == 'Filtered Signal']}

if __name__ == '__main__':
    app.run_server(debug=True)
