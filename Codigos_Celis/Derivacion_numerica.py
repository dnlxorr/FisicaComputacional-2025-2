import numpy as np
import matplotlib.pyplot as plt
import math

def f1(x):
    return math.sin(x)

def derivada(fMd,fmd,dx):
    return (fMd-fmd)/(2*dx)

def seg_derivada(fMd,fmd,f,dx):
    return (fMd-2*f+fmd)/(dx*dx)

dx=1e-1

N= 100

l_malla = N*dx


x_val = []

for i in range(N):
    x_val.append(dx*i)

y_val = []

for x in x_val:
    y_val.append(f1(x))

f_prim = []
f_prim.append(0)
for i in range(1,N-1):
    f_prim.append(derivada(y_val[i+1],y_val[i-1],dx))

f_prim.append(0)

f_dprim = []
f_dprim.append(0)
for i in range(1,N-1):
    f_dprim.append(seg_derivada(y_val[i+1],y_val[i-1],y_val[i],dx))

f_dprim.append(0)



plt.plot(x_val,y_val,color='blue')
plt.plot(x_val,f_prim,color='red')
plt.plot(x_val,f_dprim,color='green')
plt.show()

