import numpy as np
import math

def trapecios(f, a, b, n):
    h = (b - a) / n
    x = np.linspace(a, b, n+1)
    y = f(x)
    return h * (y[0] + 2*np.sum(y[1:-1]) + y[-1]) / 2

def romberg(f, a, b, tol=1e-10, max_nivel=10):
    R = np.zeros((max_nivel, max_nivel))
    for i in range(max_nivel):
        n = 2**i
        R[i, 0] = trapecios(f, a, b, n)
        for j in range(1, i+1):
            R[i, j] = (4**j * R[i, j-1] - R[i-1, j-1]) / (4**j - 1)
        if i > 0 and abs(R[i, i] - R[i-1, i-1]) < tol:
            return R[i, i]
    return R[max_nivel-1, max_nivel-1]

# Ejemplo físico: periodo del péndulo no lineal
L = 1.0
g = 9.81
theta0 = math.radians(60)
k = math.sin(theta0/2)

def integrando(phi):
    return 1.0 / np.sqrt(1.0 - (k**2) * (np.sin(phi)**2))

I = romberg(integrando, 0, math.pi/2)
T = 4 * math.sqrt(L/g) * I

print(f"Periodo del péndulo (θ_0=60°): T = {T:.12f} s")

