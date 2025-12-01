import numpy as np
import math


# FUNCIONES


def funcion1(t):
    #Función polinómica
    return t**3 + 5*t + 3

def funcion2(t):
    #Función exponencial escalada
    return 2 * math.pi * np.exp(t)



# Interpolación de Lagrange (versión local 4 nodos)


def lagrange_local(x_eval, x_nodes, y_nodes):
    resultado = 0.0
    m = len(x_nodes)

    for j in range(m):
        L = 1.0
        for k in range(m):
            if k != j:
                L *= (x_eval - x_nodes[k]) / (x_nodes[j] - x_nodes[k])
        resultado += y_nodes[j] * L
    return resultado



# Cálculo de spline cúbico natural


def obtener_coef_spline(x, y):
    n = len(x)
    h = np.diff(x)

    M = np.zeros((n, n))
    rhs = np.zeros(n)

    # Condiciones naturales
    M[0, 0] = 1
    M[-1, -1] = 1

    for i in range(1, n-1):
        M[i, i-1] = h[i-1]
        M[i, i]   = 2 * (h[i-1] + h[i])
        M[i, i+1] = h[i]

        rhs[i] = 3 * ((y[i+1] - y[i]) / h[i] -
                      (y[i] - y[i-1]) / h[i-1])

    c_vals = np.linalg.solve(M, rhs)

    a_vals = y[:-1]
    b_vals = (y[1:] - y[:-1]) / h - h * (2*c_vals[:-1] + c_vals[1:]) / 3
    d_vals = (c_vals[1:] - c_vals[:-1]) / (3*h)

    return a_vals, b_vals, c_vals[:-1], d_vals


def evaluar_spline(xq, x_tab, a, b, c, d):
    idx = np.searchsorted(x_tab, xq) - 1

    if idx < 0:
        idx = 0
    if idx >= len(a):
        idx = len(a) - 1

    dx = xq - x_tab[idx]
    return a[idx] + b[idx]*dx + c[idx]*dx**2 + d[idx]*dx**3



# PROCESADOR GENERAL DE INTERPOLACIONES


def ejecutar_interpolaciones(fun, label, n):
    print("\n===========================================")
    print(f"   EVALUANDO {label} con {n} NODOS")
    print("===========================================")

    L = 10
    paso = L / n

    xs = np.array([k*paso for k in range(n)])
    ys = fun(xs)

    x_obj = 2
    y_real = fun(x_obj)

    # ------- Lineal -------
    pos = 0
    for k in range(len(xs)-1):
        if xs[k] <= x_obj:
            pos = k
        else:
            break

    y_lin = ys[pos] + (ys[pos+1] - ys[pos]) * (x_obj - xs[pos]) / (xs[pos+1] - xs[pos])

    # Lagrange 4 puntos
    idx = np.clip(np.arange(pos-1, pos+3), 0, len(xs)-1)

    x_sub = xs[idx]
    y_sub = ys[idx]

    y_lag = lagrange_local(x_obj, x_sub, y_sub)

    # ------- Spline cúbico -------
    a, b, c, d = obtener_coef_spline(xs, ys)
    y_spl = evaluar_spline(x_obj, xs, a, b, c, d)


    print(f"x objetivo = {x_obj}")
    print(f"Valor exacto         = {y_real}")
    print("-----------------------------------------")
    print(f"Interpolación lineal   = {y_lin}")
    print(f"Interpolación Lagrange = {y_lag}")
    print(f"Spline cúbico natural  = {y_spl}")
    print("-----------------------------------------")
    print(f"Error lineal        = {abs(y_lin - y_real)}")
    print(f"Error Lagrange      = {abs(y_lag - y_real)}")
    print(f"Error Spline        = {abs(y_spl - y_real)}\n")



# EJECUCIÓN


ejecutar_interpolaciones(funcion1, "funcion1(t) = t^3 + 5t + 3", 100)
ejecutar_interpolaciones(funcion2, "funcion2(t) = 2π e^t", 100)

ejecutar_interpolaciones(funcion1, "funcion1(t) = t^3 + 9t", 1000)
ejecutar_interpolaciones(funcion2, "funcion2(t) = 2π e^t", 1000)
