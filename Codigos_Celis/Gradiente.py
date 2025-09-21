import numpy as np
import matplotlib.pyplot as plt
import math

def f1(x,y):
    return (x*x) - (y*y)

def derivada(fMd,fmd,dx):
    return (fMd-fmd)/(2*dx)

dx=1e-1
dy=1e-1
N=100

x_val = []

for i in range(N):
    x_val.append(dx*i)

y_val = []
for i in range(N):
    y_val.append(dy*i)

z= []

for i in range(N):
    fila = []
    for j in range(N):
        fila.append(f1(x_val[j],y_val[i]))
    z.append(fila)

dzdx = []
for i in range(N):
    fila = []
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
        dzdx[i][j] = (z[i][j+1]-z[i][j-1])/(2*dx)
        dzdy[i][j] = (z[i+1][j]-z[i-1][j])/(2*dy)


coord_x = []
coord_y = []
U=[]
V=[]

for i in range(N):
    for j in range(N):
        coord_x.append(x_val[i])
        coord_y.append(y_val[i])
        U.append(dzdx[i][j])
        V.append(dzdy[i][j])

plt.figure(figsize=(8,8))
plt.contour(x_val,y_val,z,levels=15,alpha=0.6)
plt.quiver(x_val,y_val,dzdx,dzdy,color='r')
plt.show()





