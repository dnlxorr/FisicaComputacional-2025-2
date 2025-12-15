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

a=0
b=20
N=1000

resultadointegral = integrar_Trapecio(a,b,N,f)

print(f"El resultado de la integral entre los intervalos {a} y {b} es aproximadamente: {resultadointegral: .0f} ")