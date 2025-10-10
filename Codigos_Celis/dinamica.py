import numpy as np
import matplotlib.pyplot as plt
import math
import random
from Particle import IonTitanio  # usa tus clases

# -------------------
# Derivada central
# -------------------
def derivada(fMd, fmd, d):
    return (fMd - fmd) / (2 * d)

# -------------------
# Interpolación bilineal
# -------------------
def interp_bilineal(x, y, X, Y, Ex, Ey, dx, dy):
    if np.isnan(x) or np.isnan(y):
        return [0.0, 0.0]

    i = int(x / dx)
    j = int(y / dy)
    nx = len(X[0])
    ny = len(Y)

    if i < 0: i = 0
    if j < 0: j = 0
    if i >= nx - 1: i = nx - 2
    if j >= ny - 1: j = ny - 2

    x0, x1 = X[0][i], X[0][i + 1]
    y0, y1 = Y[j][0], Y[j + 1][0]

    tx = (x - x0) / (x1 - x0)
    ty = (y - y0) / (y1 - y0)

    Ex_interp = ((1 - tx) * (1 - ty) * Ex[j][i] + tx * (1 - ty) * Ex[j][i + 1] +
                 (1 - tx) * ty * Ex[j + 1][i] + tx * ty * Ex[j + 1][i + 1])
    Ey_interp = ((1 - tx) * (1 - ty) * Ey[j][i] + tx * (1 - ty) * Ey[j][i + 1] +
                 (1 - tx) * ty * Ey[j + 1][i] + tx * ty * Ey[j + 1][i + 1])
    return [Ex_interp, Ey_interp]

# -------------------
# Cálculo del potencial total
# -------------------
def potencial_total(iones, X, Y):
    eps_0 = 8.854e-12
    k = 1 / (4 * math.pi * eps_0)
    Phi = np.zeros_like(X)

    for ion in iones:
        x_i, y_i = ion.getPosition()
        r = np.sqrt((X - x_i)**2 + (Y - y_i)**2)
        r[r < 1e-12] = 1e-12  # evitar singularidad
        Phi += k * ion.charge / r

    return Phi

# -------------------
# Campo eléctrico a partir del potencial
# -------------------
def campo_desde_potencial(Phi, dx, dy):
    Ex = np.zeros_like(Phi)
    Ey = np.zeros_like(Phi)

    for i in range(1, Phi.shape[1] - 1):
        for j in range(1, Phi.shape[0] - 1):
            Ex[j, i] = -derivada(Phi[j, i + 1], Phi[j, i - 1], dx)
            Ey[j, i] = -derivada(Phi[j + 1, i], Phi[j - 1, i], dy)

    return Ex, Ey

# -------------------
# Dinámica con campo derivado del potencial
# -------------------
def dinamica(iones, X, Y, dx, dy, dt, n_steps):
    for _ in range(n_steps):
        # Recalcular potencial total
        Phi = potencial_total(iones, X, Y)
        Ex, Ey = campo_desde_potencial(Phi, dx, dy)

        # Actualizar aceleraciones e integrar movimiento
        for ion in iones:
            x, y = ion.getPosition()
            ex, ey = interp_bilineal(x, y, X, Y, Ex, Ey, dx, dy)
            ax = (ion.charge / ion.mass) * ex
            ay = (ion.charge / ion.mass) * ey
            ion.setAcceleration([ax, ay])
            ion.move(dt)
            ion.applyBoundary(len(X))
    return iones

# -------------------
# Parámetros del sistema
# -------------------
Nx, Ny = 60, 60
L = 1e-6
X, Y = np.meshgrid(np.linspace(0, L, Nx), np.linspace(0, L, Ny))
dx = X[0, 1] - X[0, 0]
dy = Y[1, 0] - Y[0, 0]

# -------------------
# Inicialización de iones
# -------------------
iones = []
for _ in range(10):
    x = random.uniform(0.3e-6, 0.7e-6)
    y = random.uniform(0.3e-6, 0.7e-6)
    vx = random.uniform(-50, 50)
    vy = random.uniform(-50, 50)
    ion = IonTitanio(x, y, vx, vy)
    iones.append(ion)

# -------------------
# Simulación
# -------------------
dt = 1e-12
n_steps = 1000
iones = dinamica(iones, X, Y, dx, dy, dt, n_steps)

# -------------------
# Gráfica de trayectorias
# -------------------
plt.figure(figsize=(6, 6))
for ion in iones:
    plt.plot(ion.trayectoria_x, ion.trayectoria_y, lw=1.2)
    plt.scatter(ion.trayectoria_x[0], ion.trayectoria_y[0], color='green', marker='o')
    plt.scatter(ion.trayectoria_x[-1], ion.trayectoria_y[-1], color='red', marker='x')

plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.title('Trayectorias de iones con interacción electrostática')
plt.axis('equal')
plt.grid(True)
plt.show()
