import matplotlib.pyplot as plt
import math
import numpy as np


##metodo

def metodo_Gauss (a, b, f):
    x1 = (b - a) / 2 * (-1 / np.sqrt(3)) + (a + b) / 2
    x2 = (b - a) / 2 * (1 / np.sqrt(3)) + (a + b) / 2
    return (b - a) / 2 * (f(x1) + f(x2))
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
a_valores = [np.exp(-10),np.exp(-10), np.exp(-10), np.exp(-10), np.exp(-10)]
b_valores = [1,1,1,1,0.999]
funciones = [funcion1, funcion2, funcion3, funcion4, funcion5]
resultadosTeoricos=[1, 1, 1, math.pi/2, math.pi]


#bucle para solucionar

for i in range (len(funciones)):
    gauss = metodo_Gauss(a_valores[i], b_valores[i], funciones[i])
    print(f"\nFunción f_{i + 1}:")
    print("Cuadratura Gauss: ",gauss)
    print(f"  Valor teórico     = {resultadosTeoricos[i]}")
    print(f"  Error absoluto    = {abs(gauss - resultadosTeoricos[i]):.6f}\n")