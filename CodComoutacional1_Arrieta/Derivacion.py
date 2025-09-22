import numpy as np
import matplotlib.pyplot as plt
import math

def funcion(x):
    return math.sin(x)

def derivada(f1delta, f2delta, deltax):
    return (f1delta - f2delta) / (2 * deltax)

def segDerivada(f1delta, f2delta,f, deltax):
    return (f1delta - 2 * f + f2delta)/(deltax*deltax)

deltax = 1e-1  # = 0.01
N = 100

xValores = []
for i in range(N):
    xValores.append(i * deltax)

yValores = []
for x in xValores:
    yValores.append(funcion(x))


fprima = []
fprima.append(0)

for i in range(1, N - 1):
    fprima.append(derivada(yValores[i + 1], yValores[i - 1], deltax))

fprima.append((yValores[-1] - yValores[-2]) / deltax)


fDprim = []

fDprim.append((yValores[2] - 2*yValores[1] + yValores[0])/(deltax*deltax))
for i in range(1, N-1):
    fDprim.append(segDerivada(yValores[i+1], yValores[i-1],yValores[i],deltax))

fDprim.append((yValores[-1] - 2 * yValores[-2] + yValores[-3]) / (deltax * deltax))


x_arr = np.array(xValores)
y_arr = np.array(yValores)
fprime_arr = np.array(fprima)
f2_arr = np.array(fDprim)

print("Primeros 8 valores de f':", fprime_arr[:8])
print("Últimos 8 valores de f':", fprime_arr[-8:])

plt.plot(xValores, yValores, marker='.', linestyle='-', label='f(x)')
plt.plot(xValores, fprima, marker='.', linestyle='-', label="f'(x)")
plt.plot(xValores,fDprim, marker = '.', linestyle='-',label="f''(x)")
plt.xlabel('x')
plt.legend()
plt.grid(True)
plt.show()
