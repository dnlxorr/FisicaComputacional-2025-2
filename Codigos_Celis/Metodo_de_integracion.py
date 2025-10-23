import numpy as np
import matplotlib.pyplot as plt
import math

def integrate(a,b,N,f):
    dx = (b - a) / N
    suma = f(a)/2
    i=1

    for i in range (N-1):
        x = a+i*dx
        suma += f(x)

    suma += f(b)/2

    return suma*dx

def f (x):
    return x**3

i = integrate(0,20,1000,f)
print(f"Trapecios: {i}")

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

e = met_simpson(0,20,1000,f)
print(f"Simpson: {e}")

