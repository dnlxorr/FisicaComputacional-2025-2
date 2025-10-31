import numpy as np
import math

# Definición del metodo de Simpson 1/3

def metodo_Simpson(a, b, N, f):


    if N % 2 != 0:
        raise ValueError("N debe ser un número par para usar la regla de Simpson 1/3.")

    dx = (b - a) / N
    suma = f(a) + f(b)

    # Sumar los términos intermedios
    for i in range(1, N):
        x = a + i * dx
        if i % 2 == 0:
            suma += 2 * f(x)
        else:
            suma += 4 * f(x)

    return (dx / 3) * suma


# Ejemplo físico: Trabajo realizado por una fuerza variable F(x) = 5x²

def F(x):
    return -1  # Fuerza variable

# Intervalo de integración (desde x=0 hasta x=3)
a = 1
b = 0
N = 1000  # Debe ser par

# Cálculo del trabajo mediante el metodo de Simpson
trabajo = metodo_Simpson(a, b, N, F)

print(f"El trabajo realizado por la fuerza variable entre {a} y {b} es aproximadamente: {trabajo:.0f} J")
