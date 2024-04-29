import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, CheckButtons
from scipy.signal import cheby2, filtfilt

# Початкові значення параметрів
initial_amplitude = 1.0
initial_frequency = 1.0
initial_phase = 0.0
initial_noise_mean = 0.0
initial_noise_covariance = 0.1
initial_cutoff_frequency = 5.0

# Функція для генерації гармоніки з шумом
def harmonic_with_noise(t, amplitude, frequency, phase, noise_mean, noise_covariance, show_noise):
    signal = amplitude * np.sin(2 * np.pi * frequency * t + phase)
    if show_noise:
        noise = np.random.normal(noise_mean, np.sqrt(noise_covariance), len(t))
        return signal + noise
    else:
        return signal

# Функція для фільтрації сигналу
def filter_signal(signal, cutoff_frequency, fs=100):
    nyquist = 0.5 * fs
    cutoff = cutoff_frequency / nyquist
    b, a = cheby2(4, 40, cutoff, 'low', analog=False)
    filtered_signal = filtfilt(b, a, signal)
    return filtered_signal

# Створення основного вікна графічного інтерфейсу
fig, ax = plt.subplots()
plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.4)

t = np.arange(0.0, 10.0, 0.01)
y = harmonic_with_noise(t, initial_amplitude, initial_frequency, initial_phase, initial_noise_mean, initial_noise_covariance, True)
l, = plt.plot(t, y, lw=2, color='orange')

# Створення слайдерів
axcolor = 'lightblue'
ax_amplitude = plt.axes([0.1, 0.3, 0.65, 0.03], facecolor=axcolor)
ax_frequency = plt.axes([0.1, 0.25, 0.65, 0.03], facecolor=axcolor)
ax_phase = plt.axes([0.1, 0.2, 0.65, 0.03], facecolor=axcolor)
ax_noise_mean = plt.axes([0.1, 0.15, 0.65, 0.03], facecolor=axcolor)
ax_noise_covariance = plt.axes([0.1, 0.1, 0.65, 0.03], facecolor=axcolor)
ax_cutoff_frequency = plt.axes([0.1, 0.05, 0.65, 0.03], facecolor=axcolor)

s_amplitude = Slider(ax_amplitude, 'Amplitude', 0.1, 2.0, valinit=initial_amplitude)
s_frequency = Slider(ax_frequency, 'Frequency', 0.1, 2.0, valinit=initial_frequency)
s_phase = Slider(ax_phase, 'Phase', 0, 2*np.pi, valinit=initial_phase)
s_noise_mean = Slider(ax_noise_mean, 'Noise Mean', -1.0, 1.0, valinit=initial_noise_mean)
s_noise_covariance = Slider(ax_noise_covariance, 'Noise Covariance', 0.0, 1.0, valinit=initial_noise_covariance)
s_cutoff_frequency = Slider(ax_cutoff_frequency, 'Cutoff Frequency', 0.1, 10.0, valinit=initial_cutoff_frequency)

# Створення чекбокса
rax = plt.axes([0.8, 0.025, 0.15, 0.1], facecolor=axcolor)
check = CheckButtons(rax, ['Show Noise'], [True])

# Функція для оновлення графіку
def update(val):
    amplitude = s_amplitude.val
    frequency = s_frequency.val
    phase = s_phase.val
    noise_mean = s_noise_mean.val
    noise_covariance = s_noise_covariance.val
    cutoff_frequency = s_cutoff_frequency.val
    show_noise = check.lines[0][0].get_visible()

    y = harmonic_with_noise(t, amplitude, frequency, phase, noise_mean, noise_covariance, show_noise)
    l.set_ydata(y)

    filtered_signal = filter_signal(y, cutoff_frequency)
    l_filtered.set_ydata(filtered_signal)

    # Оновлення видимості відфільтрованого графіку в залежності від стану чекбокса
    l_filtered.set_visible(show_noise)

    fig.canvas.draw_idle()

# Призначення функції оновлення слайдерам
s_amplitude.on_changed(update)
s_frequency.on_changed(update)
s_phase.on_changed(update)
s_noise_mean.on_changed(update)
s_noise_covariance.on_changed(update)
s_cutoff_frequency.on_changed(update)

# Відфільтрована гармоніка
filtered_signal = filter_signal(y, initial_cutoff_frequency)
l_filtered, = ax.plot(t, filtered_signal, lw=2, color='purple', visible=True)

# Функція для скидання параметрів
def reset(event):
    s_amplitude.reset()
    s_frequency.reset()
    s_phase.reset()
    s_noise_mean.reset()
    s_noise_covariance.reset()
    s_cutoff_frequency.reset()

# Створення кнопки Reset
resetax = plt.axes([0.05, 0.01, 0.1, 0.02])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')
button.on_clicked(reset)

plt.show()
