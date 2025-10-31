import numpy as np
import math

# Metodo de Simpson
def metodo_Simpson(a, b, N, f):
    if N % 2 != 0:
        N += 1

    # Evitar evaluar exactamente en 0 o 1
    eps = 1e-10
    a += eps
    b -= eps

    dx = (b - a) / N
    suma = f(a)

    for i in range(1, N):
        x = a + i * dx
        if i % 2 == 0:
            suma += 2 * f(x)
        else:
            suma += 4 * f(x)

    suma += f(b)
    return suma * dx / 3


# Funciones transformadas a (0,1)
def funcion1(x):
    # Integral equivalente a e^{-x} con x = -ln(t)
    return np.exp(-(-np.log(x))) * (1/x)

def funcion2(x):
    # Evitar dominio inválido con log
    return ((np.log(1 - np.log(x))) / ((1 - np.log(x))**2)) * (1 / x)

def funcion3(x):
    return -np.log(x)

def funcion4(x):
    return (np.sin(np.log(x)) / np.log(x)) * (1/x)

def funcion5(x):
    return (1 / (1 + (np.log(x / (1 - x)))**2)) * (1 / (x * (1 - x)))


# Intervalos y funciones
a_valores = [np.exp(-10),np.exp(-10), np.exp(-10), np.exp(-10), np.exp(-10)]
b_valores = [1, 1, 1, 1, 1]
funciones = [funcion1, funcion2, funcion3, funcion4, funcion5]
resultadosTeoricos = [1, 1, 1, math.pi/2, math.pi]

# Evaluación
for i in range(len(funciones)):
    Simpson = metodo_Simpson(a_valores[i], b_valores[i], 1000, funciones[i])


    print(f"Función f_{i+1}:")
    print(f"  Método de Simpson = {Simpson}")
    print(f"  Valor teórico     = {resultadosTeoricos[i]}")
    print(f"  Error absoluto    = {abs(Simpson - resultadosTeoricos[i]):.3f}\n")
