import numpy as np
import matplotlib.pyplot as plt
import math
import random


# derivada_central: diferencia central para aproximar derivadas parciales
def derivada_central(f_plus, f_minus, dx):
    """Devuelve la derivada central (f_plus - f_minus) / (2*dx)."""
    return (f_plus - f_minus) / (2 * dx)

# pos: función auxiliar (no usada pero documentada para claridad)
def pos(x, y):
    """Retorna una tupla (x,y)."""
    return (x, y)



q_elemental = 1.602e-19
num_charges = 20
eps_0 = 8.8541e-12

k_coulomb = 1.0 / (4.0 * math.pi * eps_0)

domain_size = 1.0              # tamaño del dominio (unidad arbitraria)
n_points = 100                 # resolución de la malla (por eje)
dx = domain_size / n_points
dy = domain_size / n_points

# Coordenadas 1D para la malla
x_valores = [dx * i for i in range(n_points)]
y_valores = [dy * i for i in range(n_points)]


# Generar cargas aleatorias

charges = []  # lista de tuplas (x, y, q)
for _ in range(num_charges):
    xq = random.uniform(x_valores[0], x_valores[-1])
    yq = random.uniform(y_valores[0], y_valores[-1])
    signo = random.choice([-1, 1])
    qc = signo * q_elemental
    charges.append((xq, yq, qc))


# Construcción de mallas X,Y (listas anidadas, misma estructura que el código original)

X = []
Y = []
for j in range(n_points):
    filaX = []
    filaY = []
    for i in range(n_points):
        filaX.append(x_valores[i])
        filaY.append(y_valores[j])
    X.append(filaX)
    Y.append(filaY)


# Cálculo del potencial V(x,y)

# pot_matrix tendrá la misma estructura: pot_matrix[j][i] corresponde a (x_valores[i], y_valores[j])
pot_matrix = []
r_cutoff = 1e-2   # umbral para evitar singularidad (se omite la contribución si r <= r_cutoff)
for j in range(n_points):
    fila = []
    for i in range(n_points):
        x = X[j][i]
        y = Y[j][i]
        V_ij = 0.0
        for (xq, yq, qc) in charges:
            dx_q = x - xq
            dy_q = y - yq
            r = math.sqrt(dx_q * dx_q + dy_q * dy_q)
            if r > r_cutoff:
                V_ij += (k_coulomb * qc) / r
        fila.append(V_ij)
    pot_matrix.append(fila)


# Cálculo del campo eléctrico E = -grad(V) por diferencia central
#
Ex = [[0.0 for _ in range(n_points)] for _ in range(n_points)]
Ey = [[0.0 for _ in range(n_points)] for _ in range(n_points)]

for j in range(1, n_points - 1):
    for i in range(1, n_points - 1):
        Ex[j][i] = -derivada_central(pot_matrix[j][i + 1], pot_matrix[j][i - 1], dx)
        Ey[j][i] = -derivada_central(pot_matrix[j + 1][i], pot_matrix[j - 1][i], dy)


# Limitar magnitud del campo para visualización (evita flechas gigantes)

max_field = 1e-8
for j in range(n_points):
    for i in range(n_points):
        if abs(Ex[j][i]) > max_field:
            Ex[j][i] = math.copysign(max_field, Ex[j][i])
        if abs(Ey[j][i]) > max_field:
            Ey[j][i] = math.copysign(max_field, Ey[j][i])


# Cálculo de fuerzas netas sobre cada carga (Coulomb par a par)

forces = []  # lista de tuplas (Fx, Fy) para cada carga en 'charges'
for idx_i, (xi, yi, qi) in enumerate(charges):
    Fx, Fy = 0.0, 0.0
    for idx_j, (xj, yj, qj) in enumerate(charges):
        if idx_i != idx_j:
            dx_q = xi - xj
            dy_q = yi - yj
            r2 = dx_q * dx_q + dy_q * dy_q
            r = math.sqrt(r2)
            if r > 1e-12:
                F = k_coulomb * qi * qj / r2
                Fx += F * (dx_q / r)
                Fy += F * (dy_q / r)
    forces.append((Fx, Fy))


# Preparar vectores aplanados para quiver según el esquema solicitado
# coord_x, coord_y: listas 1D con coordenadas de cada punto de la malla (filas recorridas)
# A, B: componentes del campo eléctrico (Ex,Ey) en el mismo orden

coord_x = []
coord_y = []
A = []
B = []
for j in range(n_points):
    for i in range(n_points):
        coord_x.append(X[j][i])
        coord_y.append(Y[j][i])
        A.append(Ex[j][i])
        B.append(Ey[j][i])


# Preparar datos para graficar las cargas (ya están en 'charges')
# y para graficar las fuerzas sobre cada carga (normalizadas para mostrar dirección)

x_charges = [c[0] for c in charges]
y_charges = [c[1] for c in charges]

# Flechas de fuerza normalizadas y dimensionadas a longitud fija L (mantener esquema original)
L = 0.02
Fx_plot = []
Fy_plot = []
for fx, fy in forces:
    norm = math.hypot(fx, fy)
    if norm > 0:
        Fx_plot.append(fx / norm * L)
        Fy_plot.append(fy / norm * L)
    else:
        Fx_plot.append(0.0)
        Fy_plot.append(0.0)


# GRAFICADO (usando la plantilla que solicitaste, pero manteniendo el esquema original)
#  Contornos del potencial (pot_matrix)
# Vector field (quiver) del campo eléctrico (A,B)
# Dibujar cargas con color según signo
# Dibujar flechas en las posiciones de las cargas que indican la dirección de la fuerza neta

plt.figure(figsize=(10, 8))

# contour acepta x 1D, y 1D y z como matriz con forma (len(y), len(x))
z = pot_matrix  # alias solicitado
plt.contour(x_valores, y_valores, z, levels=15, alpha=0.8)

# quiver del campo eléctrico (flechas en escala visual)
plt.quiver(coord_x, coord_y, A, B, color='green')

# Dibujar cargas (rojo positivo, azul negativo)
for (xc, yc, qc) in charges:
    if qc > 0:
        plt.plot(xc, yc, 'ro', markersize=6)
    else:
        plt.plot(xc, yc, 'bo', markersize=6)

# Dibujar flechas de fuerza (sobre las cargas) en azul, manteniendo el estilo con pivot middle
plt.quiver(x_charges, y_charges, Fx_plot, Fy_plot,
           color='b', angles='xy', scale_units='xy', scale=1,
           width=0.008, headwidth=4, headlength=6, minlength=0,
           pivot='middle', zorder=5)

plt.xlabel("x")
plt.ylabel("y")
plt.title("Campo eléctrico (flechas) y potencial (contornos)")
plt.axis('equal')
plt.show()
