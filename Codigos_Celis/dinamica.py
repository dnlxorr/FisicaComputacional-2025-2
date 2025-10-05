import numpy as np
import matplotlib.pyplot as plt
import math
import random
from Particle import IonTitanio

# -------------------
# Derivada central
# -------------------
def derivada(fMd, fmd, dx):
    return (fMd - fmd) / (2 * dx)

# -------------------
# Interpolación bilineal (con protección de índices)
# -------------------
def interp_bilineal(x, y, X, Y, Ex, Ey, dx, dy):
    i = int(x / dx)
    j = int(y / dy)

    # Clampear índices
    nx = len(X[0])
    ny = len(Y)
    if i < 0: i = 0
    if j < 0: j = 0
    if i >= nx - 1: i = nx - 2
    if j >= ny - 1: j = ny - 2

    # Coordenadas de la celda
    x0, x1 = X[0][i], X[0][i+1]
    y0, y1 = Y[j][0], Y[j+1][0]

    # evitar división por cero
    tx = (x - x0) / (x1 - x0) if (x1 - x0) != 0 else 0.0
    ty = (y - y0) / (y1 - y0) if (y1 - y0) != 0 else 0.0

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

N = 1e-2
n_puntos = 100
dx = N / n_puntos
dy = N / n_puntos

x_val = [dx * i for i in range(n_puntos)]
y_val = [dy * i for i in range(n_puntos)]

# -------------------
# Crear varios iones de titanio
# -------------------
n_iones = 5
iones = []
for _ in range(n_iones):
    ion = IonTitanio()
    xq = random.uniform(x_val[0], x_val[-1])
    yq = random.uniform(y_val[0], y_val[-1])
    ion.setPosition([xq, yq])
    ion.setVelocity([5.0, 5.0])
    try:
        _ = ion.getAcceleration()
    except Exception:
        ion.setAcceleration([0.0, 0.0])
    iones.append(ion)

# Mallas
X = [[x_val[i] for i in range(n_puntos)] for j in range(n_puntos)]
Y = [[y_val[j] for i in range(n_puntos)] for j in range(n_puntos)]

# -------------------
# Calcular potencial
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
Ex = [[0.0 for _ in range(n_puntos)] for _ in range(n_puntos)]
Ey = [[0.0 for _ in range(n_puntos)] for _ in range(n_puntos)]

for j in range(1, n_puntos - 1):
    for i in range(1, n_puntos - 1):
        Ex[j][i] = -derivada(PHY[j][i + 1], PHY[j][i - 1], dx)
        Ey[j][i] = -derivada(PHY[j + 1][i], PHY[j - 1][i], dy)

# Limitar campo
max_field = 1e-8
for j in range(n_puntos):
    for i in range(n_puntos):
        if abs(Ex[j][i]) > max_field:
            Ex[j][i] = math.copysign(max_field, Ex[j][i])
        if abs(Ey[j][i]) > max_field:
            Ey[j][i] = math.copysign(max_field, Ey[j][i])

# -------------------
# Dinámica Leapfrog guardando posición, velocidad y aceleración
# -------------------
def dinamica(iones, X, Y, Ex, Ey, dx, dy, dt=1e-3, n_steps=1000):
    # Inicializar historial en cada ion
    for ion in iones:
        ion.hist_pos = []
        ion.hist_vel = []
        ion.hist_acc = []

    # Paso medio inicial de velocidades
    for ion in iones:
        vx, vy = ion.getVelocity()
        ax, ay = ion.getAcceleration()
        ion.setVelocity([vx + 0.5 * ax * dt, vy + 0.5 * ay * dt])

    for step in range(n_steps):
        # Actualizar posiciones
        for ion in iones:
            x, y = ion.getPosition()
            vx, vy = ion.getVelocity()
            x_new = x + vx * dt
            y_new = y + vy * dt
            ion.setPosition([x_new, y_new])
            if math.isnan(x_new) or math.isnan(y_new):
                x_new, y_new = 0.0, 0.0
            x_new=max(0.0,min(N,x_new))
            y_new=max(0.0,min(N,y_new))

            ion.setPosition([x_new, y_new])

        # Calcular aceleraciones
        for ion in iones:
            xi, yi = ion.getPosition()
            Epx, Epy = interp_bilineal(xi, yi, X, Y, Ex, Ey, dx, dy)
            Fx, Fy = ion.charge * Epx, ion.charge * Epy
            ax, ay = Fx / ion.mass, Fy / ion.mass
            ion.setAcceleration([ax, ay])

        # Actualizar velocidades
        for ion in iones:
            vx, vy = ion.getVelocity()
            ax, ay = ion.getAcceleration()
            ion.setVelocity([vx + ax * dt, vy + ay * dt])

        # Guardar datos
        for ion in iones:
            ion.hist_pos.append(tuple(ion.getPosition()))
            ion.hist_vel.append(tuple(ion.getVelocity()))
            ion.hist_acc.append(tuple(ion.getAcceleration()))

    return iones

# -------------------
# Ejecutar dinámica
# -------------------
dt = 1e-3
n_steps = 1000
iones = dinamica(iones, X, Y, Ex, Ey, dx, dy, dt=dt, n_steps=n_steps)

# -------------------
# Graficar trayectorias
# -------------------
plt.figure()
for idx, ion in enumerate(iones):
    xs = [pos[0] for pos in ion.hist_pos]
    ys = [pos[1] for pos in ion.hist_pos]
    plt.plot(xs, ys, label=f"Ion {idx}")
plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.title("Trayectorias de todos los iones")
plt.grid(True)
plt.legend()
plt.show()
