import numpy as np
import matplotlib.pyplot as plt
import math
import random

def pos(x,y):
    return (x,y)

def derivada(fMd,fmd,dx):
    return (fMd-fmd)/(2*dx)

q = 1.602e-19
n_cargas = 20
eps_0 = 8.8541e-12
k = 1.0/(4.0*math.pi*eps_0)  # constante de Coulomb

N=1.0

n_puntos = 100
dx=N/n_puntos
dy=N/n_puntos

x_val = [dx*i for i in range(n_puntos)]
y_val = [dy*i for i in range(n_puntos)]

# generar cargas aleatorias
cargas = []
for i in range(n_cargas):
    xq = random.uniform(x_val[0],x_val[-1])
    yq = random.uniform(y_val[0],y_val[-1])
    signo = random.choice([-1,1])
    qc = signo*q
    cargas.append((xq,yq,qc))

# mallas X,Y
X = []
Y = []
for j in range(n_puntos):
    filaX=[]
    filaY=[]
    for i in range(n_puntos):
        filaX.append(x_val[i])
        filaY.append(y_val[j])
    X.append(filaX)
    Y.append(filaY)

# potencial
PHY = []
cons = 1e-2
for j in range(n_puntos):
    fila =[]
    for i in range(n_puntos):
        x = X[j][i]
        y = Y[j][i]
        PHY_ij = 0.0
        for (xq,yq,qc) in cargas:
            dx_q = x-xq
            dy_q = y-yq
            r = math.sqrt(dx_q*dx_q + dy_q*dy_q)
            if r > cons:
                PHY_ij += (k * qc) / r
        fila.append(PHY_ij)
    PHY.append(fila)

# campo E
Ex = [[0.0 for j in range(n_puntos)] for i in range(n_puntos)]
Ey = [[0.0 for j in range(n_puntos)] for i in range(n_puntos)]

for j in range(1,n_puntos-1):
    for i in range(1,n_puntos-1):
        Ex[j][i] = -derivada(PHY[j][i+1],PHY[j][i-1],dx)
        Ey[j][i] = -derivada(PHY[j+1][i],PHY[j-1][i],dy)

# limitar campo
max_field = 1e-8
for j in range(n_puntos):
    for i in range(n_puntos):
        if abs(Ex[j][i])>max_field:
            Ex[j][i]=math.copysign(max_field,Ex[j][i])
        if abs(Ey[j][i])>max_field:
            Ey[j][i]=math.copysign(max_field,Ey[j][i])

# Calcular fuerza sobre cada carga
fuerzas = []
for i,(xi,yi,qi) in enumerate(cargas):
    Fx, Fy = 0.0, 0.0
    for j,(xj,yj,qj) in enumerate(cargas):
        if i != j:
            dx_q = xi - xj
            dy_q = yi - yj
            r2 = dx_q*dx_q + dy_q*dy_q
            r = math.sqrt(r2)
            if r>1e-12:  # evitar división por cero
                F = k*qi*qj/r2
                Fx += F*dx_q/r
                Fy += F*dy_q/r
    fuerzas.append((Fx,Fy))


plt.figure(figsize=(8,8))
plt.contour(X,Y,PHY,levels=15,alpha=0.6)
plt.quiver(X,Y,Ex,Ey,color='r',alpha=0.5)


x_cargas = [c[0] for c in cargas]
y_cargas = [c[1] for c in cargas]

L = 0.02
Fx_plot = []
Fy_plot = []
for fx,fy in fuerzas:
    norm = math.hypot(fx,fy)
    if norm>0:
        Fx_plot.append(fx/norm * L)
        Fy_plot.append(fy/norm * L)
    else:
        Fx_plot.append(0.0)
        Fy_plot.append(0.0)

plt.quiver(x_cargas, y_cargas, Fx_plot, Fy_plot,
           color='b', angles='xy', scale_units='xy', scale=1,
           width=0.008, headwidth=4, headlength=6, minlength=0,
           pivot='middle', zorder=5)

plt.title('Potencial, campo eléctrico y fuerza sobre cada carga')
plt.axis('equal')
plt.show()
