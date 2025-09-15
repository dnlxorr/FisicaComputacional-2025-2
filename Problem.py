import random
import numpy as np
import matplotlib.pyplot as plt

# --------------------------
# Entradas
# --------------------------
tipo = int(input("Digite 1 para enteros, 2 para flotantes: "))
a = float(input("Ingrese el límite inferior a: "))
b = float(input("Ingrese el límite superior b: "))
n = int(input("Ingrese el tamaño de cada array n: "))
m = int(input("Ingrese el número de arrays m: "))

# --------------------------
# Generar m arrays de tamaño n con números aleatorios entre [a, b]
# --------------------------
arrays = []
for _ in range(m):
    if tipo == 1:  # enteros aleatorios
        arr = [random.randint(int(a), int(b)) for _ in range(n)]
    elif tipo == 2:  # flotantes aleatorios
        arr = [random.uniform(a, b) for _ in range(n)]
    else:
        print("⚠️ Opción no válida, se generarán flotantes por defecto.")
        arr = [random.uniform(a, b) for _ in range(n)]

    arr.sort()  # orden ascendente
    arrays.append(arr)

# --------------------------
# Retomar arrays de mínimos, máximos y promedios
# --------------------------
minimos = [min(arr) for arr in arrays]
maximos = [max(arr) for arr in arrays]
promedios = [round(np.mean(arr), 2) for arr in arrays]

# --------------------------
# Mostrar resultados
# --------------------------
print("\nArrays generados (ordenados):")
for i, arr in enumerate(arrays, 1):
    print(f"Array {i}: {arr}")

print("\nArray con valores mínimos:", minimos)
print("Array con valores máximos:", maximos)
print("Array con promedios:", promedios)

# --------------------------
# Graficar mínimos y máximos
# --------------------------
plt.figure(figsize=(10, 6))
plt.plot(minimos, 'bo--', linewidth=2, markersize=8, label="Mínimos")
plt.plot(maximos, 'ro--', linewidth=2, markersize=8, label="Máximos")

plt.title("Valores mínimos y máximos de cada array")
plt.xlabel("Número de array")
plt.ylabel("Valor")
plt.legend()
plt.grid(True)
plt.show()