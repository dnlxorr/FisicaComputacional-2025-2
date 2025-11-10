import matplotlib.pyplot as plt
import math
import numpy as np


#solucion Centro de masas metodo de trapecios

def metodo_Trapecio(a,b,N,f):


    dx = (b - a) / N
    suma = f(a) / 2
    for i in range(1, N):
        x = a + i * dx
        suma += f(x)
    suma += f(b) / 2
    return suma * dx

#solucion centro de masa metodo de simpson

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

#solcion metodo de cuadraturas
def metodo_Gauss (a, b, f):
    x1 = (b - a) / 2 * (-1 / np.sqrt(3)) + (a + b) / 2
    x2 = (b - a) / 2 * (1 / np.sqrt(3)) + (a + b) / 2
    return (b - a) / 2 * (f(x1) + f(x2))

def funcion1(x):
    return 3*x**2

def funcion2(x):
    return 3*x

a_val=[0,0]
b_val=[2,2]
funciones= [funcion1,funcion2]
resultado1=[0,0]
resultado2=[0,0]
resultado3=[0,0]



for i in range(len(funciones)):
    trapecio = metodo_Trapecio(a_val[i],b_val[i], 100, funciones[i])
    Simpson = metodo_Simpson(a_val[i], b_val[i], 1000, funciones[i])
    gauss = metodo_Gauss(a_val[i], b_val[i], funciones[i])

    resultado1[i]=trapecio
    resultado2[i]=Simpson
    resultado3[i]=gauss


    print(f"Función f_{i + 1}:")
    print(f"  Método de Trapecio = {trapecio}")
    print(f" Método de Simpson = {Simpson}")
    print (f" Método de cuadraturas = {gauss}")

centroMasa1=resultado1[0]/resultado1[1]
centroMasa2=resultado2[0]/resultado2[1]
centroMasa3=resultado3[0]/resultado3[1]

print(f"Obteniendo el centro de masa por método de trapecios es: {centroMasa1 :.3f}" )
print(f"resultado metodo de Simpson= {centroMasa2 : .3f} ")
print(f"resultado metodo de Gauss = {centroMasa3 : .3f}")




