import matplotlib.pyplot as plt
import math
import numpy as np

def funcion(x,y):
    return (x*x)-(y*y)

def derivada(f1delta, f2delta, deltax):
    return (f1delta - f2delta) / (2 * deltax)

deltax= 1e-1
deltay=1e-1

N = 100

x_valores = []

for i in range (N):
    x_valores.append(i * deltax)

y_valores=[]
for i in range (N):
    y_valores.append( i * deltay)


z= []

for i in range (N):
    fila=[]
    for j in range(N):
        fila.append(funcion(x_valores[j],y_valores[i]))
    z.append(fila)


dzdx = []

for i in range (N):
    fila=[]
    for j in range(N):
        fila.append(0.0)

    dzdx.append(fila)

dzdy = []
for i in range(N):
    fila = []
    for j in range(N):
        fila.append(0.0)
    dzdy.append(fila)


for i in range(1,N-1):
    for j in range(1,N-1):
        dzdx[i][j] = (z[i][j+1]-z[i][j-1])/(2*deltax)
        dzdy[i][j] = (z[i+1][j]-z[i-1][j])/(2*deltay)

coord_x = []
coord_y = []
A=[]
B=[]

for i in range(N):
    for j in range(N):

        coord_x.append(x_valores[j])
        coord_y.append(y_valores[i])
        A.append(dzdx[i][j])
        B.append(dzdy[i][j])

plt.figure(figsize=(10,8))
plt.contour(x_valores,y_valores,z,levels=14,alpha=0.4)
plt.quiver(coord_x, coord_y, A, B, color='grey')
plt.show()