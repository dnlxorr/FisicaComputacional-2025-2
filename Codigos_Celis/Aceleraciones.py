import numpy as np
import matplotlib.pyplot as plt
import math
import random
from Particle import IonTitanio

def derivada(fMd, fmd, dx):
    return (fMd - fmd) / (2 * dx)

# ----------------------
# Interpolación bilineal
# ----------------------
def interp_bilineal(x, y, X, Y, Ex, Ey, dx, dy):
    i = int(x / dx)
    j = int(y / dy)

    # Chequeo de límites
    if i < 0 or j < 0 or i >= len(X[0]) - 1 or j >= len(Y) - 1:
        return 0.0, 0.0

    # Coordenadas de la celda
    x0, x1 = X[0][i], X[0][i+1]
    y0, y1 = Y[j][0], Y[j+1][0]

    tx = (x - x0) / (x1 - x0)
    ty = (y - y0) / (y1 - y0)

    # Ex interpolado
    Ex_interp = ((1-tx)*(1-ty)*Ex[j][i]     + tx*(1-ty)*Ex[j][i+1] +
                 (1-tx)*ty*Ex[j+1][i]       + tx*ty*Ex[j+1][i+1])

    # Ey interpolado
    Ey_interp = ((1-tx)*(1-ty)*Ey[j][i]     + tx*(1-ty)*Ey[j][i+1] +
                 (1-tx)*ty*Ey[j+1][i]       + tx*ty*Ey[j+1][i+1])

    return Ex_interp, Ey_interp

# -------------------
# Constantes
# -------------------
eps_0 = 8.8541e-12
k = 1.0 / (4.0 * math.pi * eps_0)

N = 1.0
n_puntos = 100
dx = N / n_puntos
dy = N / n_puntos

x_val = [dx * i for i in range(n_puntos)]
y_val = [dy * i for i in range(n_puntos)]

# -------------------
# Crear varios iones de titanio en posiciones aleatorias
# -------------------
n_iones = 5
iones = []
for _ in range(n_iones):
    ion = IonTitanio()
    xq = random.uniform(x_val[0], x_val[-1])
    yq = random.uniform(y_val[0], y_val[-1])
    ion.setPosition([xq, yq])
    ion.setVelocity([0.0, 0.0])
    iones.append(ion)

# Mallas
X = [[x_val[i] for i in range(n_puntos)] for j in range(n_puntos)]
Y = [[y_val[j] for i in range(n_puntos)] for j in range(n_puntos)]

# -------------------
# Calcular potencial con los iones
# -------------------
PHY = []
cons = 1e-2
for j in range(n_puntos):
    fila = []
    for i in range(n_puntos):
        x = X[j][i]
        y = Y[j][i]
        PHY_ij = 0.0
        for ion in iones:
            xq, yq = ion.getPosition()
            qc = ion.charge
            dx_q = x - xq
            dy_q = y - yq
            r = math.sqrt(dx_q * dx_q + dy_q * dy_q)
            if r > cons:
                PHY_ij += (k * qc) / r
        fila.append(PHY_ij)
    PHY.append(fila)

# -------------------
# Calcular campo
# -------------------
Ex = [[0.0 for j in range(n_puntos)] for i in range(n_puntos)]
Ey = [[0.0 for j in range(n_puntos)] for i in range(n_puntos)]

for j in range(1, n_puntos - 1):
    for i in range(1, n_puntos - 1):
        Ex[j][i] = -derivada(PHY[j][i + 1], PHY[j][i - 1], dx)
        Ey[j][i] = -derivada(PHY[j + 1][i], PHY[j - 1][i], dy)

# Limitar campo (para evitar infinitos)
max_field = 1e-8
for j in range(n_puntos):
    for i in range(n_puntos):
        if abs(Ex[j][i]) > max_field:
            Ex[j][i] = math.copysign(max_field, Ex[j][i])
        if abs(Ey[j][i]) > max_field:
            Ey[j][i] = math.copysign(max_field, Ey[j][i])

# -------------------
# Calcular fuerzas y aceleraciones con interpolación
# -------------------
pos_iones = []
Fuerzas = []
Aceleraciones = []

for ion in iones:
    xi, yi = ion.getPosition()
    Epx, Epy = interp_bilineal(xi, yi, X, Y, Ex, Ey, dx, dy)
    Fx, Fy = ion.charge * Epx, ion.charge * Epy
    ax, ay = Fx / ion.mass, Fy / ion.mass
    ion.setAcceleration([ax, ay])

    pos_iones.append((xi, yi))
    Fuerzas.append((Fx, Fy))
    Aceleraciones.append((ax, ay))

# -------------------
# Gráfico 1: Potencial + Campo
# -------------------
X_arr = np.array(X)
Y_arr = np.array(Y)
Ex_arr = np.array(Ex)
Ey_arr = np.array(Ey)

plt.figure(figsize=(7, 7))
plt.contourf(X, Y, PHY, levels=30, cmap="RdYlBu")
plt.colorbar(label="Potencial (V)")
plt.quiver(X_arr[::5, ::5], Y_arr[::5, ::5], Ex_arr[::5, ::5], Ey_arr[::5, ::5],
           color='k', alpha=0.6, scale=5e-7)
plt.title('Potencial (colormap) y Campo eléctrico (negro)')
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.axis("equal")
plt.show()

# -------------------
# Gráfico 2: Aceleraciones de los iones
# -------------------
plt.figure(figsize=(7,7))
plt.contourf(X, Y, PHY, levels=30, cmap="RdYlBu")
plt.colorbar(label="Potencial (V)")

mags_acc = [math.sqrt(ax**2 + ay**2) for (ax,ay) in Aceleraciones]
max_acc = max(mags_acc) if max(mags_acc) > 0 else 1

for (x, y), (ax, ay), mag in zip(pos_iones, Aceleraciones, mags_acc):
    plt.scatter(x, y, color="g", s=50)
    plt.quiver(x, y, (ax/max_acc)*0.3, (ay/max_acc)*0.3,
               color="b", alpha=0.7, scale=1, angles="xy", scale_units="xy")
plt.title("Aceleraciones normalizadas de los iones de Titanio")
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.axis("equal")
plt.show()

# -------------------
# Gráfico 3: Fuerzas sobre los iones
# -------------------
mags_F = [math.sqrt(Fx**2 + Fy**2) for (Fx,Fy) in Fuerzas]
max_F = max(mags_F) if max(mags_F) > 0 else 1

plt.figure(figsize=(8,8))
plt.contourf(X_arr, Y_arr, PHY, levels=30, cmap='RdYlBu')
plt.colorbar(label="Potencial (V)")

for (x, y), (Fx, Fy), mag in zip(pos_iones, Fuerzas, mags_F):
    plt.scatter(x, y, color='green', s=50)
    plt.quiver(x, y, (Fx/max_F)*0.3, (Fy/max_F)*0.3,
               color=plt.cm.viridis(mag/max_F), angles='xy', scale_units='xy', scale=1)

plt.title("Fuerzas (normalizadas) sobre los iones de Titanio")
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.axis("equal")
plt.show()

# -------------------
# Mostrar resultados
# -------------------
for idx, ion in enumerate(iones):
    print(f"Ion {idx+1}: Pos={ion.getPosition()}  Fuerza={Fuerzas[idx]}  Aceleración={ion.getAcceleration()}")







