import numpy as np
import matplotlib.pyplot as plt
import math
import random
from Particle import IonTitanio

def calcular_derivada(valor_mas, valor_menos, paso):


    return (valor_mas - valor_menos) / (2 * paso)

# -------------------
# Constantes físicas y de simulación
# -------------------
perm_vacio = 8.8541e-12  # Permitividad del vacío (F/m)
const_coulomb = 1.0 / (4.0 * math.pi * perm_vacio)  # Constante de Coulomb

tamano_dominio = 1.0  # Tamaño del dominio cuadrado (1 m)
num_puntos = 50  # Número de puntos en cada eje
delta_x = tamano_dominio / num_puntos  # Espaciado en x
delta_y = tamano_dominio / num_puntos  # Espaciado en y

coordenadas_x = [delta_x * i for i in range(num_puntos)]  # Coordenadas x
coordenadas_y = [delta_y * i for i in range(num_puntos)]  # Coordenadas y

# -------------------
# Generar iones de titanio en posiciones aleatorias
# -------------------
num_iones = 5  # Cantidad de iones
lista_iones = []
for _ in range(num_iones):
    ion = IonTitanio()
    pos_x = random.uniform(coordenadas_x[0], coordenadas_x[-1])  # Posición x aleatoria
    pos_y = random.uniform(coordenadas_y[0], coordenadas_y[-1])  # Posición y aleatoria
    ion.setPosition([pos_x, pos_y])  # Establecer posición
    ion.setVelocity([0.0, 0.0])  # Velocidad inicial cero
    lista_iones.append(ion)

# Crear mallas de coordenadas
malla_x = [[coordenadas_x[i] for i in range(num_puntos)] for j in range(num_puntos)]
malla_y = [[coordenadas_y[j] for i in range(num_puntos)] for j in range(num_puntos)]

# -------------------
# Calcular el potencial eléctrico debido a los iones
# -------------------
potencial = []
dist_min = 1e-2  # Distancia mínima para evitar singularidades
for j in range(num_puntos):
    fila_potencial = []
    for i in range(num_puntos):
        x = malla_x[j][i]
        y = malla_y[j][i]
        potencial_ij = 0.0
        for ion in lista_iones:
            x_ion, y_ion = ion.getPosition()
            carga_ion = ion.charge
            delta_x_ion = x - x_ion
            delta_y_ion = y - y_ion
            distancia = math.sqrt(delta_x_ion**2 + delta_y_ion**2)
            if distancia > dist_min:
                potencial_ij += (const_coulomb * carga_ion) / distancia
        fila_potencial.append(potencial_ij)
    potencial.append(fila_potencial)

# -------------------
# Calcular el campo eléctrico
# -------------------
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

# -------------------
# Calcular fuerzas y aceleraciones de los iones
# -------------------
posiciones_iones = []
fuerzas_iones = []
aceleraciones_iones = []

for ion in lista_iones:
    x_ion, y_ion = ion.getPosition()
    indice_x = int(x_ion / delta_x)
    indice_y = int(y_ion / delta_y)
    if 0 <= indice_x < num_puntos and 0 <= indice_y < num_puntos:
        campo_px, campo_py = campo_x[indice_y][indice_x], campo_y[indice_y][indice_x]
        fuerza_x, fuerza_y = ion.charge * campo_px, ion.charge * campo_py
        accel_x, accel_y = fuerza_x / ion.mass, fuerza_y / ion.mass
        ion.setAcceleration([accel_x, accel_y])
    else:
        fuerza_x, fuerza_y = 0.0, 0.0
        accel_x, accel_y = 0.0, 0.0
        ion.setAcceleration([0.0, 0.0])

    posiciones_iones.append((x_ion, y_ion))
    fuerzas_iones.append((fuerza_x, fuerza_y))
    aceleraciones_iones.append((accel_x, accel_y))

# -------------------
# Gráfico 1: Potencial y campo eléctrico
# -------------------
malla_x_arr = np.array(malla_x)
malla_y_arr = np.array(malla_y)
campo_x_arr = np.array(campo_x)
campo_y_arr = np.array(campo_y)

plt.figure(figsize=(7, 7))
plt.contourf(malla_x, malla_y, potencial, levels=30, cmap="RdYlBu")
plt.colorbar(label="Potencial eléctrico (V)")
plt.quiver(malla_x_arr[::5, ::5], malla_y_arr[::5, ::5], campo_x_arr[::5, ::5], campo_y_arr[::5, ::5],
           color='k', alpha=0.6, scale=5e-7)
plt.title('Potencial eléctrico (mapa de color) y Campo eléctrico (flechas negras)')
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.axis("equal")
plt.show()

# -------------------
# Gráfico 2: Aceleraciones de los iones
# -------------------
plt.figure(figsize=(7, 7))
plt.contourf(malla_x, malla_y, potencial, levels=30, cmap="RdYlBu")
plt.colorbar(label="Potencial eléctrico (V)")

magnitudes_acel = [math.sqrt(ax**2 + ay**2) for (ax, ay) in aceleraciones_iones]
max_acel = max(magnitudes_acel) if max(magnitudes_acel) > 0 else 1

for (x, y), (ax, ay), mag in zip(posiciones_iones, aceleraciones_iones, magnitudes_acel):
    plt.scatter(x, y, color="g", s=50)
    plt.quiver(x, y, (ax/max_acel)*0.3, (ay/max_acel)*0.3,
               color="b", alpha=0.7, scale=1, angles="xy", scale_units="xy")
plt.title("Aceleraciones normalizadas de los iones de titanio")
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.axis("equal")
plt.show()

# -------------------
# Gráfico 3: Fuerzas sobre los iones
# -------------------
magnitudes_fuerza = [math.sqrt(fx**2 + fy**2) for (fx, fy) in fuerzas_iones]
max_fuerza = max(magnitudes_fuerza) if max(magnitudes_fuerza) > 0 else 1

plt.figure(figsize=(8, 8))
plt.contourf(malla_x_arr, malla_y_arr, potencial, levels=30, cmap='RdYlBu')
plt.colorbar(label="Potencial eléctrico (V)")

for (x, y), (fx, fy), mag in zip(posiciones_iones, fuerzas_iones, magnitudes_fuerza):
    plt.scatter(x, y, color='green', s=50)
    plt.quiver(x, y, (fx/max_fuerza)*0.3, (fy/max_fuerza)*0.3,
               color=plt.cm.viridis(mag/max_fuerza), angles='xy', scale_units='xy', scale=1)

plt.title("Fuerzas normalizadas sobre los iones de titanio")
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.axis("equal")
plt.show()

# -------------------
# Mostrar resultados finales
# -------------------
for idx, ion in enumerate(lista_iones):
    print(f"Ion {idx+1}: Posición={ion.getPosition()}  Fuerza={fuerzas_iones[idx]}  Aceleración={ion.getAcceleration()}")