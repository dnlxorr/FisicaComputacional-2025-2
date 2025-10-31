import math
import numpy as np
import matplotlib.pyplot as plt

# Métodos numéricos
def met_trapecios(a, b, N, f):
    dx = (b - a) / N
    suma = f(a) / 2
    for i in range(1, N):
        x = a + i * dx
        suma += f(x)
    suma += f(b) / 2
    return suma * dx

def met_simpson(a, b, N, f):
    if N % 2 != 0:
        N += 1
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

def met_cuadraturas_Gauss(a, b, f):
    x1 = (b - a) / 2 * (-1 / np.sqrt(3)) + (a + b) / 2
    x2 = (b - a) / 2 * (1 / np.sqrt(3)) + (a + b) / 2
    return (b - a) / 2 * (f(x1) + f(x2))

# Funciones
def f_1(x):
    return 1

def f_2(x):
    return np.exp(-(np.log(x))**2) / x

def f_3(x):
    return 2 * np.cos(np.log(x)) * np.exp(-(np.log(x))**2) / x

def f_4(x):
    return 2 / (1 + (np.log(x))**2) * (1/x)

def f_5(x):
    return (np.log(1/x)) * x

# Intervalos
a_vals = [0, np.exp(-10), np.exp(-10), np.exp(-10), np.exp(-10)]
b_vals = [0.1, 1, 1, 1, 1]
funcs = [f_1, f_2, f_3, f_4, f_5]
valores_teoricos = [0.1, np.sqrt(np.pi)/2, 1.38039, math.pi, 0.25]

# Bucle general
for i in range(len(funcs)):
    trap = met_trapecios(a_vals[i], b_vals[i], 1000, funcs[i])
    simp = met_simpson(a_vals[i], b_vals[i], 1000, funcs[i])
    gauss = met_cuadraturas_Gauss(a_vals[i], b_vals[i], funcs[i])

    print(f"\nFunción f_{i+1}:")
    print("Método de Trapecios: ",trap)
    print("Método de Simpson: ",simp)
    print("Cuadratura Gauss: ",gauss)
    print("Valor teórico: ",valores_teoricos[i])
    print(f"Error Trapecios: {abs(trap - valores_teoricos[i]):.6f}")
    print(f"Error Simpson: {abs(simp - valores_teoricos[i]):.6f}")
    print(f"Error Gauss: {abs(gauss - valores_teoricos[i]):.6f}")


