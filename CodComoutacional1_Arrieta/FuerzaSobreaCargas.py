import numpy as np
import matplotlib.pyplot as plt
import math
import random
import particle

def obtener_posicion(x, y):
    """
    Devuelve una tupla con las coordenadas (x, y).
    """
    return (x, y)

def calcular_derivada(valor_mas, valor_menos, paso):
    """
    Calcula la derivada numérica usando el método de diferencia central.
    """
    return (valor_mas - valor_menos) / (2 * paso)

# -------------------
# Constantes físicas y de simulación
# -------------------

carga_elemental = 1.602e-19  # Carga elemental (C)
num_cargas = 20
perm_vacio = 8.8541e-12
const_coulomb = 1.0 / (4.0 * math.pi * perm_vacio)

tamano_dominio = 1.0

num_puntos = 100
delta_x = tamano_dominio / num_puntos
delta_y = tamano_dominio / num_puntos

coordenadas_x = [delta_x * i for i in range(num_puntos)]  # Coordenadas x
coordenadas_y = [delta_y * i for i in range(num_puntos)]  # Coordenadas y


# Generar cargas puntuales aleatorias

lista_cargas = []
for _ in range(num_cargas):
    pos_x = random.uniform(coordenadas_x[0], coordenadas_x[-1])  # Posición x aleatoria
    pos_y = random.uniform(coordenadas_y[0], coordenadas_y[-1])  # Posición y aleatoria
    signo = random.choice([-1, 1])  # Signo aleatorio (+ o -)
    valor_carga = signo * carga_elemental
    lista_cargas.append((pos_x, pos_y, valor_carga))


# Crear mallas de coordenadas

malla_x = []
malla_y = []
for j in range(num_puntos):
    fila_x = []
    fila_y = []
    for i in range(num_puntos):
        fila_x.append(coordenadas_x[i])
        fila_y.append(coordenadas_y[j])
    malla_x.append(fila_x)
    malla_y.append(fila_y)


# Calcular el potencial eléctrico

potencial = []
dist_min = 1e-2  # Distancia mínima para evitar singularidades
for j in range(num_puntos):
    fila_potencial = []
    for i in range(num_puntos):
        x = malla_x[j][i]
        y = malla_y[j][i]
        potencial_ij = 0.0
        for (x_ion, y_ion, carga_ion) in lista_cargas:
            delta_x_ion = x - x_ion
            delta_y_ion = y - y_ion
            distancia = math.sqrt(delta_x_ion**2 + delta_y_ion**2)
            if distancia > dist_min:
                potencial_ij += (const_coulomb * carga_ion) / distancia
        fila_potencial.append(potencial_ij)
    potencial.append(fila_potencial)


# Calcular el campo eléctrico

campo_x = [[0.0 for j in range(num_puntos)] for i in range(num_puntos)]
campo_y = [[0.0 for j in range(num_puntos)] for i in range(num_puntos)]

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


# Calcular la fuerza sobre cada carga

fuerzas_iones = []
for i, (x_i, y_i, carga_i) in enumerate(lista_cargas):
    fuerza_x, fuerza_y = 0.0, 0.0
    for j, (x_j, y_j, carga_j) in enumerate(lista_cargas):
        if i != j:  # Evitar autointeracción
            delta_x_ion = x_i - x_j
            delta_y_ion = y_i - y_j
            distancia_cuadrada = delta_x_ion**2 + delta_y_ion**2
            distancia = math.sqrt(distancia_cuadrada)
            if distancia > 1e-12:  # Evitar división por cero
                fuerza = const_coulomb * carga_i * carga_j / distancia_cuadrada
                fuerza_x += fuerza * delta_x_ion / distancia
                fuerza_y += fuerza * delta_y_ion / distancia
    fuerzas_iones.append((fuerza_x, fuerza_y))


# Visualización: Potencial, campo eléctrico y fuerzas

plt.figure(figsize=(8, 8))
plt.contour(malla_x, malla_y, potencial, levels=15, alpha=0.6)
plt.quiver(malla_x, malla_y, campo_x, campo_y, color='g', alpha=0.5)

posiciones_x_cargas = [carga[0] for carga in lista_cargas]
posiciones_y_cargas = [carga[1] for carga in lista_cargas]

longitud_flecha = 0.02
fuerzas_x_grafico = []
fuerzas_y_grafico = []
for fx, fy in fuerzas_iones:
    norma = math.hypot(fx, fy)
    if norma > 0:
        fuerzas_x_grafico.append(fx / norma * longitud_flecha)
        fuerzas_y_grafico.append(fy / norma * longitud_flecha)
    else:
        fuerzas_x_grafico.append(0.0)
        fuerzas_y_grafico.append(0.0)

plt.quiver(posiciones_x_cargas, posiciones_y_cargas, fuerzas_x_grafico, fuerzas_y_grafico,
           color='b', angles='xy', scale_units='xy', scale=1,
           width=0.008, headwidth=4, headlength=6, minlength=0,
           pivot='middle', zorder=5)

plt.title('Potencial eléctrico, campo eléctrico y fuerzas sobre las cargas')
plt.axis('equal')
plt.show()