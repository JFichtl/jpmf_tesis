import numpy as np
from scipy.ndimage import convolve
from matplotlib import pyplot as plt

class kernel: 
    def __init__(self, name, a, b, n, points=1000): 
        self.name = name 
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
    
        A = b+a
        B = b-a
        f   = 2*np.pi/B
        phi = np.pi*A/B
        self.x = f*self.t - phi
        if self.name == 'fejer': 
            self.values  = (np.sin((n+1)* self.x *0.5) / np.sin(0.5* self.x ) )**2 / (B*(1+n))
        elif self.name == 'dirichlet':
            self.values =  (1/B) * np.sin((n+0.5)* self.x)/np.sin(0.5* self.x )
        elif self.name == 'poisson':
            # add a regularization constant to avoid division by zero,
            # define a sequence n that goes to 1 as n grows
            r = 1 - (1/(n+1)) + 1e-6
            self.values =  (1/B) * (1-pow(r,2))/(1-2*r*np.cos(self.x) + pow(r,2))
        elif self.name == 'heat':
            t = 1/(n+1)
            self.values = (1/pow(4*np.pi*t, 0.5)) * np.exp(-pow(self.x/pow(4*t, 0.5) , 2))


    def convolve(self, array, show = False, null=False, mode = 'circular'):
        #t_conv  = np.linspace(2*self.a, 2*self.b, 2*len(array) -1)
        if mode == 'circular':
            t_conv = self.t
            convolution = self.delta * convolve(self.values, array, mode='wrap')
        elif mode == 'linear':
            t_conv  = np.linspace(2*self.a, 2*self.b, 2*len(array) -1)
            convolution = self.delta * np.convolve(self.values, array)
        if show:
            plt.plot( self.t, array, label = 'Initial f.')
            if not null:
                plt.plot( t_conv,  convolution, label = f"{self.n}th order approx.")
                plt.legend()
                plt.show()
        else:
            return t_conv, convolution
    def show(self):
        plt.plot(self.t, self.values)
        plt.show()

def approximate_by_kernel(array, a, b, start=1, end=10, skip=1, ker_name='dirichlet'):
    mode = 'circular'
    if ker_name in ['dirichlet', 'fejer', 'poisson', 'heat']:
        if ker_name == 'heat':
            mode = 'linear'
        kernel(ker_name, a, b, 0).convolve(array, show=True, null = True, mode=mode)
        for i in range(start, end+skip, skip):
           Ker = kernel(ker_name, a, b, i)
           x, y = Ker.convolve(array, mode=mode)
           plt.plot(x, y, label = f"k = {i}", alpha = 1 - 0.5*(i/end) )
        plt.legend()
        plt.show()
    else:
        raise("Invalid kernel")
