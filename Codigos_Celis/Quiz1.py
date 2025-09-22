import numpy as np
import matplotlib.pyplot as plt
import math
import random

def pos(x,y):
    return (x,y)

def derivada(fMd,fmd,dx):
    return (fMd-fmd)/(2*dx)

q = 1.602e-19
n_cargas = 2
eps_0 = 8.8541e-12
k = 1.0/(4.0*math.pi*eps_0)

N=1.0

n_puntos = 100
dx=N/n_puntos
dy=N/n_puntos

x_val = []

for i in range(n_puntos):
    x_val.append(dx*i)

y_val = []
for i in range(n_puntos):
    y_val.append(dy*i)

cargas = []
for i in range(n_cargas):
    xq = random.uniform(x_val[0],x_val[-1])
    yq = random.uniform(y_val[0],y_val[-1])
    signo = random.choice([-1,1])
    qc = signo*q
    cargas.append((xq,yq,qc))


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

Ex = [] #dpdx
for i in range(n_puntos):
    fila = []
    for j in range(n_puntos):
        fila.append(0.0)
    Ex.append(fila)

Ey = [] #dpdy
for i in range(n_puntos):
    fila = []
    for j in range(n_puntos):
        fila.append(0.0)
    Ey.append(fila)


for j in range(1,n_puntos-1):
    for i in range(1,n_puntos-1):
        Ex[j][i] = -derivada(PHY[j][i+1],PHY[j][i-1],dx)
        Ey[j][i] = -derivada(PHY[j+1][i],PHY[j-1][i],dy)



print(cargas)
'''
for j in range (n_puntos):
    for i in range(n_puntos):
        norm = math.sqrt(Ex[j][i]*Ex[j][i]+Ey[j][i]*Ey[j][i])
        if norm != 0:
            Ex[j][i] /= norm
            Ey[j][i] /= norm
'''
max_field = 1e-8  # valor máximo permitido
for j in range(n_puntos):
    for i in range(n_puntos):
        if abs(Ex[j][i])>max_field:
            Ex[j][i]=math.copysign(max_field,Ex[j][i])
        if abs(Ey[j][i])>max_field:
            Ey[j][i]=math.copysign(max_field,Ey[j][i])

plt.figure(figsize=(8,8))
plt.contour(X,Y,PHY,levels=15,alpha=0.6)
plt.quiver(X,Y,Ex,Ey,color='r')

plt.title('Potencial y campo electrico')
plt.show()
