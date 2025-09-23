import random
import math
import numpy as np
import matplotlib.pyplot as plt

#  Parámetros
L = 1       # tamaño del cuadrado
N = 10           # número de cargas
nx =80        # malla en x
ny =80        # malla en y
k = 1.0
eps = 1e-1
delta = 1e-2
random.seed(123)

#  distribuir cargas
charges = []
for _ in range(N):
    x_i = random.random() * L
    y_i = random.random() * L
    sign = random.choice([-1, 1])
    q_i = sign * (0.5 + random.random()*0.5)
    charges.append((x_i, y_i, q_i))

#  malla y potencial
x_valores = [i * (L / (nx - 1)) for i in range(nx)]
y_valores = [j * (L / (ny - 1)) for j in range(ny)]

phi = [[0.0 for _ in range(nx)] for _ in range(ny)]

for i in range(ny):
    y = y_valores[i]
    for j in range(nx):
        x = x_valores[j]
        v = 0.0
        for (xc, yc, qc) in charges:

            dx = x - xc
            dy = y - yc
            r = math.sqrt(dx*dx + dy*dy) + eps
            if r < delta:
                pass
            else:
                v += k * qc / r
        phi[i][j] = v

# gradiente
dx = x_valores[1] - x_valores[0]
dy = y_valores[1] - y_valores[0]

dzdx = [[0.0 for _ in range(nx)] for _ in range(ny)]
dzdy = [[0.0 for _ in range(nx)] for _ in range(ny)]

for i in range(ny):
    for j in range(nx):
        if 0 < j < nx - 1:
            dphidx = (phi[i][j+1] - phi[i][j-1]) / (2 * dx)
        elif j == 0:
            dphidx = (phi[i][j+1] - phi[i][j]) / dx
        else:
            dphidx = (phi[i][j] - phi[i][j-1]) / dx

        if 0 < i < ny - 1:
            dphidy = (phi[i+1][j] - phi[i-1][j]) / (2 * dy)
        elif i == 0:
            dphidy = (phi[i+1][j] - phi[i][j]) / dy
        else:
            dphidy = (phi[i][j] - phi[i-1][j]) / dy

        dzdx[i][j] = -dphidx
        dzdy[i][j] = -dphidy


coord_x = []
coord_y = []
A = []
B = []
for i in range(ny):
    for j in range(nx):
        coord_x.append(x_valores[j])
        coord_y.append(y_valores[i])
        A.append(dzdx[i][j])
        B.append(dzdy[i][j])

# Convertir phi en array para contour
z = np.array(phi)

plt.figure(figsize=(10,8))
plt.contour(x_valores, y_valores, z, levels=15, alpha=0.8)
plt.quiver(coord_x, coord_y, A, B, color='grey')

# Dibujar cargas
for (xc, yc, qc) in charges:
    if qc > 0:
        plt.plot(xc, yc, 'ro', markersize=6)
    else:
        plt.plot(xc, yc, 'bo', markersize=6)

plt.xlabel("x")
plt.ylabel("y")
plt.title("Campo eléctrico (flechas) y potencial (contornos)")
plt.show()
