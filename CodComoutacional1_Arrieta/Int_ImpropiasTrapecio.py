import matplotlib.pyplot as plt
import math
import numpy as np

## metodo

def metodo_Trapecio(a,b,N,f):


    dx = (b - a) / N
    suma = f(a) / 2
    for i in range(1, N):
        x = a + i * dx
        suma += f(x)
    suma += f(b) / 2
    return suma * dx

##Funciones
def funcion1 (x):
    return np.exp(-x)

def funcion2 (x):
    return ((np.log(1 - np.log(x))) / ((1 - np.log(x))**2)) * (1/x)

def funcion3 (x):
    return -np.log(x)

def funcion4(x):
    return (np.sin(np.log(x))/np.log(x))*1/x

def funcion5(x):
    return (1 / (1 + (np.log(x/1-x))**2)) * (1 / x * (1-x))

##intervalos
a_valores = [np.exp(-10),np.exp(-10), np.exp(-10), np.exp(-10), np.exp(-2)]
b_valores = [1,1,1,0.99999,0.99999]
funciones = [funcion1, funcion2, funcion3, funcion4, funcion5]
resultadosTeoricos=[1,1,1,math.pi/2,math.pi]

for i in range(len(funciones)):
    trapecio = metodo_Trapecio(a_valores[i], b_valores[i], 1000, funciones[i])
    print(f"Función f_{i+1}:")
    print(f"  Método de Trapecio = {trapecio}")
    print(f"  Valor teórico     = {resultadosTeoricos[i]}")
    print(f"  Error absoluto    = {abs(trapecio - resultadosTeoricos[i]):.6f}\n")