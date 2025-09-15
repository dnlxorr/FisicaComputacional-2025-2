import numpy as np
import matplotlib.pyplot as plt

# Entrada de datos
n = int(input("Ingrese el tamaño de cada array (n): "))
m = int(input("Ingrese el número de arrays (m): "))

print("Tipo de números:")
print("1. Enteros")
print("2. Decimales")
tipo = int(input("Seleccione el tipo (1 o 2): "))

min_val = float(input("Ingrese el valor mínimo del intervalo: "))
max_val = float(input("Ingrese el valor máximo del intervalo: "))

# Generación de la matriz de NumPy
if tipo == 1:
    # Para enteros, se usa np.random.randint
    arrays_np = np.random.randint(int(min_val), int(max_val) + 1, size=(m, n))
else:
    # Para decimales, se usa np.random.uniform
    arrays_np = np.random.uniform(min_val, max_val, size=(m, n))

# Ordenar cada fila de la matriz
for i in range(m):
    arrays_np[i].sort()

# Visualización de la matriz en un mapa de calor
plt.imshow(arrays_np, cmap='coolwarm', interpolation='nearest')
plt.colorbar(label='Valor')
plt.title('Mapa de Calor de los Arrays Ordenados')
plt.xlabel('Índice dentro del Array')
plt.ylabel('Número de Array')
plt.show()

# --- Funciones de cálculo ---

def prom_arrays(arr_np):
    # Calcula el promedio de cada fila (eje 1)
    return np.mean(arr_np, axis=1)

def val_may(arr_np):
    # El valor máximo de un array ordenado es el último elemento
    # Se obtienen los valores de la última columna
    return arr_np[:, -1]

# --- Cálculo y visualización de promedios y máximos ---
lista_mayores = val_may(arrays_np)
promedio = prom_arrays(arrays_np)

indices = range(1, len(lista_mayores) + 1)

plt.figure(figsize=(8, 5))
plt.plot(indices, lista_mayores, 'o-r', label='Máximos')
plt.plot(indices, promedio, 's-b', label='Promedios')

plt.xlabel('Número de array')
plt.ylabel('Valor')
plt.title('Valores máximos y promedios de cada array')
plt.legend()
plt.grid(True)
plt.show()