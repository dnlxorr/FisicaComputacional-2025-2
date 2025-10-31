
import numpy as np
import matplotlib.pyplot as plt
import math
import random
from particle import IonTitanio


# Calcular derivada central

def calcular_derivada(valor_mas, valor_menos, paso):

    #Calcula la derivada numérica usando el método de diferencia central.

    return (valor_mas - valor_menos) / (2 * paso)


# Interpolación bilineal con protección de índices

def interpolacion_bilineal(x, y, malla_x, malla_y, campo_x, campo_y, delta_x, delta_y):

    #Realiza interpolación bilineal para obtener los valores del campo eléctrico en (x, y).

    indice_x = int(x / delta_x)
    indice_y = int(y / delta_y)


    # Asegurar que los índices estén dentro de los límites

    num_x = len(malla_x[0])
    num_y = len(malla_y)
    if indice_x < 0: indice_x = 0
    if indice_y < 0: indice_y = 0
    if indice_x >= num_x - 1: indice_x = num_x - 2
    if indice_y >= num_y - 1: indice_y = num_y - 2


    # Coordenadas de la celda de la malla

    x0, x1 = malla_x[0][indice_x], malla_x[0][indice_x + 1]
    y0, y1 = malla_y[indice_y][0], malla_y[indice_y + 1][0]


    # Calcular factores de interpolación, evitando división por cero

    tx = (x - x0) / (x1 - x0) if (x1 - x0) != 0 else 0.0
    ty = (y - y0) / (y1 - y0) if (y1 - y0) != 0 else 0.0


    # Interpolación para campo_x

    campo_x_interp = ((1 - tx) * (1 - ty) * campo_x[indice_y][indice_x] +
                      tx * (1 - ty) * campo_x[indice_y][indice_x + 1] +
                      (1 - tx) * ty * campo_x[indice_y + 1][indice_x] +
                      tx * ty * campo_x[indice_y + 1][indice_x + 1])

    # Interpolación para campo_y

    campo_y_interp = ((1 - tx) * (1 - ty) * campo_y[indice_y][indice_x] +
                      tx * (1 - ty) * campo_y[indice_y][indice_x + 1] +
                      (1 - tx) * ty * campo_y[indice_y + 1][indice_x] +
                      tx * ty * campo_y[indice_y + 1][indice_x + 1])

    return campo_x_interp, campo_y_interp


# Constantes físicas y de simulación

perm_vacio = 8.8541e-12  # Permitividad del vacío (F/m)
const_coulomb = 1.0 / (4.0 * math.pi * perm_vacio)  # Constante de Coulomb

tamano_dominio = 1e-2  # Tamaño del dominio cuadrado (m)
num_puntos = 100  # Número de puntos en cada eje
delta_x = tamano_dominio / num_puntos  # Espaciado en x
delta_y = tamano_dominio / num_puntos  # Espaciado en y

coordenadas_x = [delta_x * i for i in range(num_puntos)]  # Coordenadas x
coordenadas_y = [delta_y * i for i in range(num_puntos)]  # Coordenadas y

# Generar iones de titanio en posiciones aleatorias


num_iones = 5  # Cantidad de iones
lista_iones = []
for _ in range(num_iones):
    ion = IonTitanio()
    pos_x = random.uniform(coordenadas_x[0], coordenadas_x[-1])  # Posición x aleatoria
    pos_y = random.uniform(coordenadas_y[0], coordenadas_y[-1])  # Posición y aleatoria
    ion.establecer_posicion([pos_x, pos_y])
    ion.establecer_velocidad([5.0, 5.0])  # Velocidad inicial fija
    try:
        _ = ion.obtener_aceleracion()
    except Exception:
        ion.establecer_aceleracion([0.0, 0.0])
    lista_iones.append(ion)


# Crear mallas de coordenadas

malla_x = [[coordenadas_x[i] for i in range(num_puntos)] for j in range(num_puntos)]
malla_y = [[coordenadas_y[j] for i in range(num_puntos)] for j in range(num_puntos)]


# Calcular el potencial eléctrico

potencial = []
dist_min = 1e-2  # Distancia mínima para evitar singularidades
for j in range(num_puntos):
    fila_potencial = []
    for i in range(num_puntos):
        x = malla_x[j][i]
        y = malla_y[j][i]
        potencial_ij = 0.0
        for ion in lista_iones:
            x_ion, y_ion = ion.obtener_posicion()
            carga_ion = ion.carga
            delta_x_ion = x - x_ion
            delta_y_ion = y - y_ion
            distancia = math.sqrt(delta_x_ion**2 + delta_y_ion**2)
            if distancia > dist_min:
                potencial_ij += (const_coulomb * carga_ion) / distancia
        fila_potencial.append(potencial_ij)
    potencial.append(fila_potencial)


# Calcular el campo eléctrico

campo_x = [[0.0 for _ in range(num_puntos)] for _ in range(num_puntos)]
campo_y = [[0.0 for _ in range(num_puntos)] for _ in range(num_puntos)]

for j in range(1, num_puntos - 1):
    for i in range(1, num_puntos - 1):
        campo_x[j][i] = -calcular_derivada(potencial[j][i + 1], potencial[j][i - 1], delta_x)
        campo_y[j][i] = -calcular_derivada(potencial[j + 1][i], potencial[j - 1][i], delta_y)


# Limitar el campo eléctrico para evitar valores extremos

limite_campo = 1e-8
for j in range(num_puntos):
    for i in range(num_puntos):
        if abs(campo_x[j][i]) > limite_campo:
            campo_x[j][i] = math.copysign(limite_campo, campo_x[j][i])
        if abs(campo_y[j][i]) > limite_campo:
            campo_y[j][i] = math.copysign(limite_campo, campo_y[j][i])


# Simulación dinámica usando el método Leapfrog

def simulacion_dinamica(iones, malla_x, malla_y, campo_x, campo_y, delta_x, delta_y, paso_tiempo=1e-3, num_pasos=1000):

    # Inicializar historiales para cada ion
    for ion in iones:
        ion.historial_posiciones = []
        ion.historial_velocidades = []
        ion.historial_aceleraciones = []

    # Actualizar velocidades con medio paso inicial
    for ion in iones:
        vel_x, vel_y = ion.obtener_velocidad()
        acel_x, acel_y = ion.obtener_aceleracion()
        ion.establecer_velocidad([vel_x + 0.5 * acel_x * paso_tiempo, vel_y + 0.5 * acel_y * paso_tiempo])

    for paso in range(num_pasos):
        # Actualizar posiciones
        for ion in iones:
            x, y = ion.obtener_posicion()
            vel_x, vel_y = ion.obtener_velocidad()
            x_nuevo = x + vel_x * paso_tiempo
            y_nuevo = y + vel_y * paso_tiempo
            if math.isnan(x_nuevo) or math.isnan(y_nuevo):
                x_nuevo, y_nuevo = 0.0, 0.0
            x_nuevo = max(0.0, min(tamano_dominio, x_nuevo))
            y_nuevo = max(0.0, min(tamano_dominio, y_nuevo))

            ion.establecer_posicion([x_nuevo, y_nuevo])

        # Calcular aceleraciones usando interpolación bilineal
        for ion in iones:
            x_i, y_i = ion.obtener_posicion()
            campo_px, campo_py = interpolacion_bilineal(x_i, y_i, malla_x, malla_y, campo_x, campo_y, delta_x, delta_y)
            fuerza_x, fuerza_y = ion.carga * campo_px, ion.carga * campo_py
            acel_x, acel_y = fuerza_x / ion.masa, fuerza_y / ion.masa
            ion.establecer_aceleracion([acel_x, acel_y])

        # Actualizar velocidades
        for ion in iones:
            vel_x, vel_y = ion.obtener_velocidad()
            acel_x, acel_y = ion.obtener_aceleracion()
            ion.establecer_velocidad([vel_x + acel_x * paso_tiempo, vel_y + acel_y * paso_tiempo])

        # Guardar datos en los historiales
        for ion in iones:
            ion.historial_posiciones.append(tuple(ion.obtener_posicion()))
            ion.historial_velocidades.append(tuple(ion.obtener_velocidad()))
            ion.historial_aceleraciones.append(tuple(ion.obtener_aceleracion()))

    return iones


# Ejecutar la simulación dinámica

paso_tiempo = 1e-3  # Intervalo de tiempo (s)
num_pasos = 1000  # Número de pasos de simulación
lista_iones = simulacion_dinamica(lista_iones, malla_x, malla_y, campo_x, campo_y, delta_x, delta_y, paso_tiempo=paso_tiempo, num_pasos=num_pasos)


# Visualizar trayectorias de los iones

plt.figure()
for idx, ion in enumerate(lista_iones):
    posiciones_x = [pos[0] for pos in ion.historial_posiciones]
    posiciones_y = [pos[1] for pos in ion.historial_posiciones]
    plt.plot(posiciones_x, posiciones_y, label=f"Ion {idx}")
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.title("Trayectorias de los iones de titanio")
plt.grid(True)
plt.legend()
plt.show()












