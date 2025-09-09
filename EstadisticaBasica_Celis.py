import numpy as np
import random

n = int(input('ingrese el tamaño de cada array: '))
m = int(input('ingrese el número de arrays que desea: '))

arrays = []
while True:
    print('Seleccione el valor que desea para los números dentro de los arrays')
    print('1. int (entero)')
    print('2. float (decimal)')
    print('3. Salir ')
    op = int(input())

    if op ==3:
        break

    if op == 1:
        print('Ingrese desde donde comienza el intervalo:')
        a= int(input())
        print('Ingrese donde termina el intervalo:')
        b= int(input())

        for i in range(m):
            x_i = np.random.randint(a,b, size=n)
            arrays.append(x_i)
        break

    elif op == 2:
        print('Ingrese desde donde comienza el intervalo:')
        a = float(input())
        print('Ingrese donde termina el intervalo:')
        b = float(input())

        for i in range(m):
            x_i = np.random.uniform(a, b, size=n)
            arrays.append(x_i)
        break
    else:
        print('Error: ingrese una opción válida')

def ordenar_arrays (arrays):
    for i in arrays:
        i.sort()

def prom_arrays (arrays):
    promedio = []
    for i in arrays:
        x= i.mean()
        promedio.append(x)
    return promedio

def val_may (arrays):
    ordenar_arrays(arrays)
    mayores = []
    for i in arrays:
        x = i[len(i)-1]
        mayores.append(x)
    return mayores

ordenar_arrays(arrays)

lista_mayores = val_may(arrays)
promedio = prom_arrays(arrays)

print(lista_mayores)
print(promedio)