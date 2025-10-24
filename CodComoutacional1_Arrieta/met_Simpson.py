
import numpy as np
import math

# Definimos el metodo de simpson 1/3

def metodo_Simpson(a, b, N, f):

    if N % 2 != 0:
        raise ValueError("N debe ser par para aplicar la regla de Simpson 1/3.")

    dx = (b - a) / N
    suma = f(a) + f(b)

    # Bucle para sumar los términos interiores
    for i in range(1, N):
        x = a + i * dx
        if i % 2 == 0:
            suma += 2 * f(x)
        else:
            suma += 4 * f(x)

    return (dx / 3) * suma

# Ejemplo físico: cálculo del trabajo realizado por una fuerza variable F(x) = 2x
# entre x = 0 y x = 10
def F(x):
    return 2 * x **2

# Llamamos al metodo
resultado = metodo_Simpson(0, 10, 1000, F)

print("El trabajo aproximado (integral de F(x)) es:", resultado)
