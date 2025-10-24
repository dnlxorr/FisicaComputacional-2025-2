import numpy as np
import matplotlib.pyplot as plt
import math

def integrar_Trapecio(a,b,N,f):
    dx = (b - a) / N
    suma = f(a)/2
    i=1

    for i in range (N-1):
        x = a+i*dx
        suma += f(x)

    suma += f(b)/2

    return suma*dx

def f (x):
    return x

i = integrate(0,20,1000,f)
print(i)