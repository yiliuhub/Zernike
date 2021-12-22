import torch
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt 

class Zernike:
    def __init__(self, order=6, reflocxy=[], deltaxy=[]):
        self.order = order 
        self.reflocxy = reflocxy 
        self.deltaxy = deltaxy
        print(reflocxy, deltaxy)


    def monoterm(self, powerx, powery, x, y):
        return np.power(x, powerx) * np.power(y, powery)


order = 2
z = Zernike(order, reflocxy=1, deltaxy=1)

print(z.monoterm(2, 2, 2, 3))
