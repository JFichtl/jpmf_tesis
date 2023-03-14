import numpy as np
from scipy.ndimage import convolve
from matplotlib import pyplot as plt

trad = {'dirichlet':'Dirichlet', 'fejer':'Fejér', 'poisson':'Poisson', 'heat':'calor'}

def hamming_dist(w1, w2):
    A = set(w1)
    B = set(w2)
    n = len(A.intersection(B))
    m = len(A.symmetric_difference(A))
    k = len(A.union(B))
    return n/k

class kernel: 
    def __init__(self, name, a, b, n, points=1000): 
        self.name = name.lower()
        self.n    = n
        self.a    = a
        self.b    = b
        if n < 0:
            raise Exception("n can't be less than zero.")
        if a>=b:
            raise Exception("a has to be strictly less than 1.")
        if (b-a)<= 0.01:
            raise Exception("a and b are less than the minimum precision.")
        self.t    = np.linspace(a, b, points)
        # delta is needed to approximate the integral
        self.delta = self.t[1] - self.t[0]
    
        if self.name in trad.keys():
           self.values = self.ker_generator(self.n) 
        else:
            comparison = [hamming_dist(self.name, w2) for w2 in trad.keys() ]
            m = max(comparison)
            likely = list(trad.keys())[comparison.index(m)]
            print(f"I'm guessing you meant {likely}.")
            self.name = likely
            self.values = self.ker_generator(self.n)
            
    def ker_generator(self, n):
        A = self.b + self.a
        B = self.b - self.a
        f   = 2*np.pi/B
        phi = np.pi*A/B
        x = f * self.t - phi
        if self.name == 'fejer':
          return (np.sin((n + 1)* x *0.5) / np.sin(0.5* x ) )**2 / (B*(1+n))
        elif self.name == 'dirichlet':
           return  (1/B) * np.sin((n+0.5)* x)/np.sin(0.5* x )
        elif self.name == 'poisson':
            # add a regularization constant to avoid division by zero,
            # define a sequence n that goes to 1 as n grows
            r = 1 - (1/(n+1)) + 1e-6
            return  (1/B) * (1-pow(r,2))/(1-2*r*np.cos(x) + pow(r,2))
        elif self.name == 'heat':
            t = 1/(n+1)
            return (2*np.pi/B)*(1/pow(4*np.pi*t, 0.5)) * np.exp(-pow(x/pow(4*t, 0.5) , 2))


    def convolve(self, array, show = False, null=False, mode = 'circular'):
        #t_conv  = np.linspace(2*self.a, 2*self.b, 2*len(array) -1)
        if mode == 'circular':
            t_conv = self.t
            convolution = self.delta * convolve(self.values, array, mode='wrap')
        elif mode == 'linear':
            t_conv  = np.linspace(2*self.a, 2*self.b, 2*len(array) -1)
            convolution = self.delta * np.convolve(self.values, array)
        if show:
            plt.plot( self.t, array, label = 'Función inicial f.')
            if not null:
                plt.plot( t_conv,  convolution, label = f"Aproximación de orden {self.n}.")
                plt.legend()
                plt.title(f"Convolución del núcleo de {trad[self.name]} con f.")
                plt.show()
        else:
            return t_conv, convolution
            
    def approximate(self, array, start=1, end=10, skip=1):
        plt.clf()
        mode = 'circular'
        if self.name == 'heat':
            mode = 'linear'
        
        self.convolve(array, show=True, null = True, mode=mode)
        for i in range(start, end+skip, skip):
            self.values = self.ker_generator(i)
            x, y = self.convolve(array, mode=mode)
            plt.plot(x, y, label = f"k = {i}", alpha = 1 - 0.5*(i/end) )
        self.values = self.ker_generator(self.n)
        plt.title(f"Aproximación de f usando el núcleo de {trad[self.name]}")
        plt.legend()
        plt.show()
        
    def show(self):
        plt.plot(self.t, self.values)
        plt.show()

def showcase_kernel(a=-np.pi, b=np.pi, start=1, end=10, skip=1, ker_name = 'dirichlet'):
    if ker_name in ['dirichlet', 'fejer', 'poisson', 'heat']:
        for i in range(start, end+skip, skip):
            k = kernel(ker_name, a, b, i)
            x, y = k.t, k.values
            plt.plot(x, y, label = f'n={i}')
        plt.title(f'Núcleo de {trad[ker_name]}')
        plt.legend()
        plt.show()

