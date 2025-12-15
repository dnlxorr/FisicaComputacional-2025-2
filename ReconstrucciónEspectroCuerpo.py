#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
reconstruccion_spectral_simplificado.py

Versión simplificada para estudiantes (pregrado).
Genera espectros sintéticos (Ley de Planck), muestrea, añade ruido,
reconstruye con CubicSpline, estima lambda_max, calcula integral parcial,
grafica y guarda resultados.

Requisitos: numpy, scipy, matplotlib, pandas
"""

import os
from math import isfinite
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq, curve_fit
import pandas as pd

# ----------------------- CONFIGURACIÓN (CAMBIA AQUÍ) -----------------------
CONFIG = {
    "experiment_name": "exp_simple",   # prefijo de archivos de salida
    "T": 5800,                         # temperatura en Kelvin
    "lmin_um": 0.1,                    # lambda mínima en micrómetros
    "lmax_um": 3.0,                    # lambda máxima en micrómetros
    "N": 25,                           # número de puntos muestreados
    "pattern": "uniform",              # 'uniform' | 'log' | 'peak_concentrated'
    "noise_type": "gaussian",          # 'gaussian' | 'poisson' | None
    "noise_alpha": 0.05,               # nivel de ruido relativo (ej. 0.05 = 5%)
    "do_planck_fit": True,             # si True, intenta ajustar A * Planck(l,T)
    "save_fig": True,                  # guarda figura en carpeta results/
    "save_csv": True,                  # guarda resumen csv en results/
    "plot_show": False                 # True para mostrar la figura (bloqueante)
}
# ---------------------------------------------------------------------------

# Constantes físicas
h = 6.62607015e-34
c = 299792458.0
kB = 1.380649e-23
WIEN_B = 2.897771955e-3  # m*K

# Directorio de resultados
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(SCRIPT_DIR, "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

# ----------------------- Funciones físicas y utilidades --------------------

def planck_lambda(l_m, T):
    """
    Ley de Planck por longitud de onda (metros).
    l_m: array o escalar (metros)
    T: temperatura (K)
    Devuelve B_lambda en SI (W sr^-1 m^-3)
    """
    l = np.array(l_m, dtype=float)
    with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
        x = (h * c) / (l * kB * T)
        denom = np.expm1(x)   # exp(x) - 1
        B = (2 * h * c**2) / (l**5 * denom)
    B = np.where(np.isfinite(B), B, 0.0)
    return B

def sample_spectrum(T, l_min_um, l_max_um, N, pattern='uniform'):
    """
    Devuelve:
      l_um: longitudes de onda en micrómetros (array length N)
      l_m: longitudes en metros
      I_true: radiancia teórica en cada punto (sin ruido)
    pattern: 'uniform', 'log', 'peak_concentrated'
    """
    if pattern == 'uniform':
        l_um = np.linspace(l_min_um, l_max_um, N)
    elif pattern == 'log':
        l_um = np.exp(np.linspace(np.log(l_min_um), np.log(l_max_um), N))
    elif pattern == 'peak_concentrated':
        # Más puntos cerca de lambda_max teórico
        lmax_theo_um = WIEN_B / T * 1e6
        n_peak = max(2, N // 2)
        n_rest = N - n_peak
        peak_range = (max(l_min_um, 0.8 * lmax_theo_um), min(l_max_um, 1.2 * lmax_theo_um))
        peak_points = np.linspace(peak_range[0], peak_range[1], n_peak)
        rest_points = np.linspace(l_min_um, l_max_um, n_rest)
        l_um = np.sort(np.concatenate([peak_points, rest_points]))
    else:
        raise ValueError("pattern debe ser 'uniform', 'log' o 'peak_concentrated'")

    l_m = l_um * 1e-6
    I_true = planck_lambda(l_m, T)
    return l_um, l_m, I_true

def add_noise(I_true, noise_type='gaussian', rel_sigma=0.05, seed=None):
    """
    Agrega ruido a I_true:
      - gaussian: ruido aditivo normal con sigma = rel_sigma * max(I_true)
      - poisson: ruido tipo conteo (aproximado)
      - None: no agrega ruido
    """
    if noise_type is None:
        return np.array(I_true, copy=True)
    rng = np.random.default_rng(seed)
    if noise_type == 'gaussian':
        sigma = rel_sigma * np.max(I_true)
        noise = rng.normal(0.0, sigma, size=I_true.shape)
        I_obs = I_true + noise
        I_obs = np.where(I_obs > 0, I_obs, 0.0)
        return I_obs
    elif noise_type == 'poisson':
        # Escalado para evitar conteos muy pequeños en intensidades SI
        scale = 1e6 / np.max(I_true)
        counts = rng.poisson(I_true * scale)
        return counts / scale
    else:
        raise ValueError("noise_type debe ser 'gaussian', 'poisson' o None")

# ----------------------- Reconstrucción y estimadores ----------------------

def reconstruct_cubic(lambda_m, I_obs):
    """Spline cúbico natural (usa lambda en metros)."""
    return CubicSpline(lambda_m, I_obs, bc_type='natural')

def find_lambda_max_from_spline(spline, lam_min, lam_max, num_intervals=1000):
    """
    Encuentra lambda_max buscando raíces de s'(lambda)=0 y eligiendo la raíz
    que sea máximo (s''<0) y mayor valor de s.
    Devuelve lambda en metros.
    """
    s1 = spline.derivative(1)
    s2 = spline.derivative(2)
    lam_grid = np.linspace(lam_min, lam_max, num_intervals)
    vals1 = s1(lam_grid)
    roots = []
    for i in range(len(lam_grid)-1):
        a, b = lam_grid[i], lam_grid[i+1]
        if np.isnan(vals1[i]) or np.isnan(vals1[i+1]):
            continue
        if vals1[i] == 0:
            roots.append(lam_grid[i])
        elif vals1[i] * vals1[i+1] < 0:
            try:
                r = brentq(lambda x: s1(x), a, b)
                roots.append(r)
            except Exception:
                pass
    # Filtrar raíces y tomar la que sea máximo local (s''<0) y mayor valor
    candidates = []
    for r in roots:
        if (r > lam_min) and (r < lam_max):
            try:
                if s2(r) < 0:
                    candidates.append((r, spline(r)))
            except Exception:
                pass
    if candidates:
        lam_max_est, _ = max(candidates, key=lambda t: t[1])
        return lam_max_est
    # fallback: máxima en malla fina
    lam_fine = np.linspace(lam_min, lam_max, max(2000, num_intervals))
    vals = spline(lam_fine)
    idx = np.nanargmax(vals)
    return lam_fine[idx]

def integrate_trapz(fun_or_spline, a, b, num=1000):
    x = np.linspace(a, b, num)
    if callable(fun_or_spline):
        y = fun_or_spline(x)
    else:
        y = fun_or_spline(x)
    mask = np.isfinite(y)
    if np.count_nonzero(mask) < 2:
        return 0.0
    return np.trapz(y[mask], x[mask])

def planck_fit(lambda_m, I_obs, T_guess=None, scale_guess=None):
    """
    Ajuste simple I_obs ≈ A * B_lambda(lambda, T). Devuelve (I_fit_func, (T_fit, A_fit))
    """
    if T_guess is None:
        idx = np.nanargmax(I_obs)
        l_guess = lambda_m[idx]
        T_guess = WIEN_B / (l_guess if l_guess > 0 else 1e-6)
    if scale_guess is None:
        scale_guess = np.max(I_obs) / (np.max(planck_lambda(lambda_m, T_guess)) + 1e-30)

    def model(l, T, A):
        return A * planck_lambda(l, T)

    try:
        popt, pcov = curve_fit(model, lambda_m, I_obs, p0=[T_guess, scale_guess], maxfev=20000)
        T_fit, A_fit = popt
        def I_fit(l): return A_fit * planck_lambda(l, T_fit)
        return I_fit, (T_fit, A_fit)
    except Exception as e:
        # si falla, devolver None
        return None, (np.nan, np.nan)

# ----------------------- Experimento simplificado --------------------------

def run_simple_experiment(cfg):
    """
    Ejecuta un experimento simple según la configuración cfg (diccionario).
    Retorna un diccionario con resultados y rutas de archivos generados.
    """
    # Generar datos (sin ruido)
    l_um, l_m, I_true = sample_spectrum(cfg["T"], cfg["lmin_um"], cfg["lmax_um"], cfg["N"], pattern=cfg["pattern"])
    # Añadir ruido (si aplica)
    I_obs = add_noise(I_true, noise_type=cfg["noise_type"], rel_sigma=cfg["noise_alpha"], seed=42)

    lam_min, lam_max = l_m.min(), l_m.max()

    # Reconstrucción spline cúbica
    spline = reconstruct_cubic(l_m, I_obs)

    # Encontrar lambda_max por spline
    lammax_spline = find_lambda_max_from_spline(spline, lam_min, lam_max)
    T_from_wien = WIEN_B / lammax_spline if lammax_spline is not None else np.nan

    # Integrales (área parcial en el rango medido)
    area_true = integrate_trapz(lambda x: planck_lambda(x, cfg["T"]), lam_min, lam_max)
    area_spline = integrate_trapz(spline, lam_min, lam_max)

    # Ajuste a Planck (opcional)
    planck_result = None
    if cfg.get("do_planck_fit", False):
        I_fit_func, params = planck_fit(l_m, I_obs)
        T_fit = params[0]
        area_planck_fit = integrate_trapz(I_fit_func, lam_min, lam_max) if I_fit_func is not None else np.nan
        planck_result = {"T_fit": T_fit, "area_fit": area_planck_fit}
    else:
        planck_result = {"T_fit": np.nan, "area_fit": np.nan}

    # --- Gráfica ---
    lam_fine = np.linspace(lam_min, lam_max, 800)
    lam_fine_um = lam_fine * 1e6

    fig, ax = plt.subplots(figsize=(8,4.5))
    ax.scatter(l_um, I_obs, label="Datos (muestreados)", s=20)
    ax.plot(l_um, I_true, label="Teórico (Planck)", linewidth=1)
    ax.plot(lam_fine_um, spline(lam_fine), label="Spline cúbico (reconstrucción)", linewidth=1)
    if planck_result["T_fit"] and not np.isnan(planck_result["T_fit"]):
        ax.plot(lam_fine_um, planck_fit(l_m, I_obs)[0](lam_fine), '--', label=f"Ajuste Planck (T≈{planck_result['T_fit']:.0f} K)")
    # marcar lambdas máximos
    ax.axvline((WIEN_B/cfg["T"]*1e6) if False else lammax_spline*1e6, linestyle='--', label="λ_max (spline)", color='k', alpha=0.6)
    ax.set_xlabel("λ (µm)")
    ax.set_ylabel("Radiancia B_λ (SI)")
    ax.set_title(f"Reconstrucción (T={cfg['T']} K, N={cfg['N']}, ruido={cfg['noise_alpha']*100:.1f}%)")
    ax.legend(fontsize='small', loc='upper right')
    fig.tight_layout()

    # Guardar figura y CSV resumen
    fig_path = None
    csv_path = None
    if cfg.get("save_fig", True):
        fig_name = f"{cfg['experiment_name']}_recon.png"
        fig_path = os.path.join(RESULTS_DIR, fig_name)
        fig.savefig(fig_path, dpi=200)
        if not cfg.get("plot_show", False):
            plt.close(fig)

    # Guardar resumen
    summary = {
        "experiment": cfg["experiment_name"],
        "T_true": cfg["T"],
        "N": cfg["N"],
        "pattern": cfg["pattern"],
        "noise_type": cfg["noise_type"],
        "noise_alpha": cfg["noise_alpha"],
        "lam_max_spline_um": lammax_spline*1e6 if isfinite(lammax_spline) else np.nan,
        "T_from_wien": T_from_wien,
        "area_true": area_true,
        "area_spline": area_spline,
        "T_fit_planck": planck_result["T_fit"],
        "area_fit_planck": planck_result["area_fit"]
    }
    if cfg.get("save_csv", True):
        csv_name = f"{cfg['experiment_name']}_summary.csv"
        csv_path = os.path.join(RESULTS_DIR, csv_name)
        df = pd.DataFrame([summary])
        if os.path.exists(csv_path):
            df_old = pd.read_csv(csv_path)
            df = pd.concat([df_old, df], ignore_index=True)
        df.to_csv(csv_path, index=False)

    return {"summary": summary, "fig_path": fig_path, "csv_path": csv_path}

# ------------------------------ EJECUCIÓN ----------------------------------
if __name__ == "__main__":
    cfg = CONFIG
    print("Ejecutando experimento simplificado con configuración:")
    for k,v in cfg.items():
        print(f"  {k}: {v}")
    result = run_simple_experiment(cfg)
    print("Resultado (resumen):")
    for k,v in result["summary"].items():
        print(f"  {k}: {v}")
    if result["fig_path"]:
        print("Figura guardada en:", result["fig_path"])
    if result["csv_path"]:
        print("CSV guardado en:", result["csv_path"])
    if cfg.get("plot_show", False):
        plt.show()
