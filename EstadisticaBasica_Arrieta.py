import random
import matplotlib.pyplot as plt

#  Entrada de datos
n = int(input("Ingrese el tamaño de cada array (n): "))
m = int(input("Ingrese el número de arrays (m): "))

print("Tipo de números:")
print("1. Enteros")
print("2. Decimales")
tipo = int(input("Seleccione el tipo (1 o 2): "))

min_val = float(input("Ingrese el valor mínimo del intervalo: "))
max_val = float(input("Ingrese el valor máximo del intervalo: "))

# Generación de arrays
arrays = []
for i in range(m):
    if tipo == 1:  # Enteros
        arr = [random.randint(int(min_val), int(max_val)) for _ in range(n)]
    else:  # Decimales
        arr = [round(random.uniform(min_val, max_val), 2) for _ in range(n)]

    arr.sort()  # ordenar ascendentemente
    arrays.append(arr)

# Cálculo de promedios
medias = []
for arr in arrays:
    media = sum(arr) / len(arr)
    medias.append(media)

#  Mostrar resultados
print("\nArrays generados y ordenados:")
for i, arr in enumerate(arrays, start=1):
    print(f"Array {i}: {arr}")

print("\nMedias de cada array:")
for i, media in enumerate(medias, start=1):
    print(f"Array {i} -> media = {media}")

#  Graficar arrays
for i, arr in enumerate(arrays, start=1):
    plt.plot(arr, marker='o', label=f'Array {i}')
plt.title('Valores de los arrays')
plt.xlabel('Índice')
plt.ylabel('Valor')
plt.legend()
plt.grid(True)
plt.show()

#  Graficar promedios
plt.bar(range(1, m + 1), medias)
plt.title('promedios de cada array')
plt.xlabel('Array')
plt.ylabel('promedio')
plt.show()
