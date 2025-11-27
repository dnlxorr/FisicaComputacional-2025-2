import math
import numpy as np

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



# Un resorte con una constante de elasticidad k = 200 N/m está en su posición
# de equilibrio. ¿Cuánto trabajo (en Joules) se requiere para estirarlo desde
# su equilibrio (x=0) hasta una longitud de x = 0.5 metros?

# La ley de Hooke es F(x) = kx.
# El trabajo es la integral de la fuerza a lo largo de la distancia.

def F (x):
    k = 200 #N/m
    return -k*x

trap = met_trapecios(0, 0.5, 1000, F)
simp = met_simpson(0, 0.5, 1000, F)
gauss = met_cuadraturas_Gauss(0, 0.5,F)

valor_teorico = -25

print(f"\nFunción F(x)=kx:")
print("Método de Trapecios: ",trap)
print("Método de Simpson: ",simp)
print("Cuadratura Gauss: ",gauss)
print("Valor teórico: ",valor_teorico)
print(f"Error Trapecios: {abs(trap - valor_teorico):.6f}")
print(f"Error Simpson: {abs(simp - valor_teorico):.6f}")
print(f"Error Gauss: {abs(gauss - valor_teorico):.6f}")

