import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from random import random as rand
precision = 0.01

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams["figure.figsize"] = (10,8)

def fourier_coefs_k(f, domain, k):
    A = 0
    B = 0
    domf = np.arange(domain[0], domain[1], precision)
    T = domain[1] - domain[0]
    omega = 2*np.pi/T
    for t in domf:
        A += f(t)*np.cos(omega*k*t)*precision
        B += f(t)*np.sin(omega*k*t)*precision
    if k == 0:
        factor = 1/T
    else:
        factor = 2/T
    A *= factor
    B *= factor
    return A, B

def fourier_approx(f, domainA, domainB, steps, plot_steps = False, plot_skip = 1):
    domf     = np.arange(domainA[0], domainA[1], precision)
    ext_domf = np.arange(domainB[0], domainB[1], precision)
    T = domainA[1] - domainA[0]
    omega = 2*np.pi/T
    f_hat = np.zeros(len(ext_domf))
    legends = ["original"]
    plt.plot(domf, f(domf))
    for k in range(0, steps + plot_skip):
        A, B = fourier_coefs_k(f, domainA, k)
        f_hat += A*np.cos(omega*k*ext_domf) + B*np.sin(omega*k*ext_domf)
        if plot_steps and (k%plot_skip == 0):
            plt.plot(ext_domf, f_hat)
            legends.append(f"k = {k}")
    if not plot_steps:
        legends.append(f"k={steps}")
        plt.plot(ext_domf, f_hat)
    plt.legend(legends)
    plt.show()

def plot_lollipop_freqs(buff_len):
    freqs     = []
    freq_vals = []
    for i in range(0, buff_len+1):
        freqs.append(rf"$w_{i}$")
        freq_vals.append(rand())
    for i in range(1, buff_len):
        freqs.append(rf"$-w_{buff_len-i}$")
        freq_vals.append(freq_vals[buff_len-i])
    plt.stem(freq_vals)
    plt.xticks(range(2*buff_len), freqs)
    plt.title("Representación Gráfica de la Salida de FFTW3")
    plt.xlabel("Frecuencia")
    plt.ylabel("Magnitud")
    plt.show()

# Aquí podemos verificar que la matriz circulante cumple:
# circulant_matrix(u,v)*v = circulant_matrix(v,u)*v
# (reemplazamos a * por np.matmul(A, B)
def circulant_matrix(u,v):
    len_u = len(u)
    len_v = len(v)
    circ_uv = np.zeros((len_u,len_v))
    for i in range(0, len_u):
        for j in range(0, len_v):
            circ_uv[i][j] = u[(len_v+i-j)%len_v]
    return circ_uv
