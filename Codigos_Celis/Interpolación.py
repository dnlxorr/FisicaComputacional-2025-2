import numpy as np
import math

# ------------------------------------------------------
# FUNCIONES
# ------------------------------------------------------

def f1(x):
    return x**3 + 9*x

def f2(x):
    return 2*math.pi*np.exp(x)

# ------------------------------------------------------
# Interpolación de Lagrange (local con 4 puntos)
# ------------------------------------------------------

def lagrange_interp(x, xs, ys):
    total = 0.0
    n = len(xs)
    for j in range(n):
        Lj = 1.0
        for k in range(n):
            if j != k:
                Lj *= (x - xs[k]) / (xs[j] - xs[k])
        total += ys[j] * Lj
    return total

# ------------------------------------------------------
# Cálculo de coeficientes del spline cúbico natural
# ------------------------------------------------------

def cubic_spline_coeffs(x, y):
    n = len(x)
    h = np.diff(x)

    A = np.zeros((n, n))
    b = np.zeros(n)

    A[0,0] = 1
    A[-1,-1] = 1

    for i in range(1, n-1):
        A[i,i-1] = h[i-1]
        A[i,i]   = 2*(h[i-1] + h[i])
        A[i,i+1] = h[i]
        b[i] = 3*((y[i+1]-y[i])/h[i] - (y[i]-y[i-1])/h[i-1])

    c = np.linalg.solve(A, b)

    a = y[:-1]
    b2 = (y[1:] - y[:-1])/h - h*(2*c[:-1] + c[1:])/3
    d = (c[1:] - c[:-1])/(3*h)

    return a, b2, c[:-1], d


def spline_eval(x, x_vals, a, b, c, d):
    # localizar intervalo
    i = np.searchsorted(x_vals, x) - 1
    if i < 0:
        i = 0
    if i >= len(a):
        i = len(a)-1

    h = x - x_vals[i]
    return a[i] + b[i]*h + c[i]*h**2 + d[i]*h**3


def procesar_interpolaciones(func, nombre, n_puntos):
    print("\n===========================================")
    print(f"  PROCESANDO {nombre} CON {n_puntos} PUNTOS")
    print("===========================================")

    N = 10
    dx = N/n_puntos
    x_vals = np.array([i*dx for i in range(n_puntos)])
    y_vals = func(x_vals)

    x0 = 2
    y_exacto = func(x0)

    # --- Lineal ---
    i = 0
    for k in range(len(x_vals)-1):
        if x_vals[k] <= x0:
            i = k
        else:
            break

    y_lin = y_vals[i] + (y_vals[i+1]-y_vals[i]) * (x0-x_vals[i])/(x_vals[i+1]-x_vals[i])

    # --- Lagrange con 4 puntos ---
    idxs = np.clip(np.arange(i-1, i+3), 0, len(x_vals)-1)
    xL = x_vals[idxs]
    yL = y_vals[idxs]

    y_lagrange = lagrange_interp(x0, xL, yL)

    # --- Spline cúbico ---
    a, b, c, d = cubic_spline_coeffs(x_vals, y_vals)
    y_spline = spline_eval(x0, x_vals, a, b, c, d)

    # --- RESULTADOS ---
    print(f"x0 = {x0}")
    print(f"Valor analítico      = {y_exacto}")
    print("-----------------------------------------")
    print(f"Interpolación lineal = {y_lin}")
    print(f"Interpolación Lagrange = {y_lagrange}")
    print(f"Spline cúbico natural = {y_spline}")
    print("-----------------------------------------")
    print(f"Error lineal        = {abs(y_lin - y_exacto)}")
    print(f"Error Lagrange      = {abs(y_lagrange - y_exacto)}")
    print(f"Error Spline        = {abs(y_spline - y_exacto)}")
    print("\n")


procesar_interpolaciones(f1, "f1(x) = x^3 + 9x", 100)
procesar_interpolaciones(f2, "f2(x) = 2π e^x", 100)

procesar_interpolaciones(f1, "f1(x) = x^3 + 9x", 1000)
procesar_interpolaciones(f2, "f2(x) = 2π e^x", 1000)

