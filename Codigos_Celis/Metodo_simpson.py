import numpy as np
import math


def met_simpson (a,b,N,f):
    dx = (b - a) / N
    suma = f(a)
    i = 1
    for i in range (N-1):
        x = a+i*dx
        if i % 2 == 0:
            suma += 2*f(x)

        if i%2 != 0:
            suma += 4*f(x)

    suma += f(b)
    return suma*(dx/3)

def f (x):
    return x**3

i = met_simpson(0,20,1000,f)
print(i)