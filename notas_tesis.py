import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from random import random as rand
precision = 0.01

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams["figure.figsize"] = (10,8)

class kernel: 
    def __init__(self, name, n, domain): 
        self.name = name 
        self.n    = n
        if n < 0:
            raise Exception("n can't be less than zero")
        self.t    = domain 
        if self.name == 'fejer': 
            self.values  = ( np.sin((n+1)*t*0.5)/ np.sin(0.5*t) )**2 / (2*np.pi*(1+n))
        elif self.name == 'dirichlet': 
            self.values  = ( np.sin((n+0.5)*t)/np.sin(0.5*t) )/(2*np.pi)
            self.signum  = np.abs(self.values) / (self.values + 1e-16)
            self.zeros   = np.linspace( domain[0],  domain[-1], 2*(2*n+1)+1)
        elif self.name == 'poisson': 
            r = 1 - 1/(n+1) 
            self.values  = (1-r**2) / (2*np.pi*(1-2*r*np.cos(t) + r**2))
        else:
            raise Exception("Invalid kernel")

# Hacemos una función específica de graficación para el núcleo de Dirichlet
# para obtener sus envolventes.
def plot_dirichlet_kernel(n = 15, skip = 3, abs_value = False, evolute = False):
    plt.clf()
    t = np.linspace(-np.pi, np.pi, 1000)
    legends = [] 
    if evolute:
        x = (np.sin(t/2))**(-1)
        plt.plot(t, x)
        plt.plot(t, -x)
        legends += [r"$y = sen(t/2)^{-1}$", r"$y = -sen(t/2)^{-1}$"]
    
    for i in range(1, n + skip, skip):
        if abs_value:
            x = np.abs(np.sin((i+0.5)*t)/np.sin(0.5*t))
        else:
            x = np.sin((i+0.5)*t)/np.sin(t/2)
        plt.plot(t, x)
        legends.append(f"k = {i}")
    #plt.title("Núcleo de Dirichlet")
    plt.legend(legends, loc = 'upper right')
    plt.show()

def plot_kernel(name, start = 0, steps = 5, skip = 3):
    plt.clf()
    legends = []
    t = np.linspace(-np.pi, np.pi, 1000)
    for i in range(start, steps*skip, skip):
        ker = kernel(name, i, t)
        x = ker.values
        plt.plot(t, x)
        legends.append(f"n = {i}")
    plt.legend(legends)
    plt.show()

def plot_kernel_convolution(domain, data, kernel_name, n = 16, skip= 3, fading_factor = 0.6, alpha_increase = 1):
    legends = [r"$f(t)$"]
    plt.plot(domain, data)
    a = domain[0]
    b = domain[-1]
    domain_length = len(domain)
    domain_conv = np.linspace(2*a, 2*b, 2*domain_length -1)
    delta_t = np.diff(domain)[0]
    alpha_index = 1
    for i in range(0, 51, 10):
        kern = kernel(kernel_name, i, domain)
        y = np.convolve(data, kern.values)*delta_t
        plt.plot(domain_conv, y, alpha = 1/pow(alpha_index, fading_factor))
        legends.append(f"$n = {i}$")
        alpha_index += alpha_increase
    plt.legend(legends)
    plt.show()

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
