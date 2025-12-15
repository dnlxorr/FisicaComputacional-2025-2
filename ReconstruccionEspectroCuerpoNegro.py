#!/usr/bin/env python3
"""
reconstruccion_spectral_ready.py

Versión corregida y lista para Windows/PyCharm.
Genera espectros sintéticos (Planck), muestrea, añade ruido,
reconstruye con CubicSpline, UnivariateSpline (smoothing) y ajuste a Planck,
estima lambda_max, temperatura vía Wien, integra el área parcial y hace bootstrap.

Salida: carpeta ./results/ con imágenes y CSV resumen.

Requisitos: numpy, scipy, matplotlib, pandas
"""

import os
import warnings
from math import isfinite

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline, UnivariateSpline
from scipy.optimize import brentq, curve_fit
import pandas as pd

# ----------------------- Constantes físicas --------------------------------
h = 6.62607015e-34
c = 299792458.0
kB = 1.380649e-23
WIEN_B = 2.897771955e-3  # m*K

# ----------------------- Utilidades de archivo ------------------------------
# Carpeta de salida (results dentro del directorio del script)
try:
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
except NameError:
    SCRIPT_DIR = os.getcwd()
RESULTS_DIR = os.path.join(SCRIPT_DIR, "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

def save_fig_safe(fig, fname):
    path = os.path.join(RESULTS_DIR, fname)
    try:
        fig.savefig(path, dpi=200)
    except Exception as e:
        print("Error guardando figura en", path, ":", e)
        raise
    finally:
        plt.close(fig)
    return path

# ----------------------- Ley de Planck -------------------------------------
def planck_lambda(l_m, T):
    """
    Ley de Planck por longitud de onda (metros).
    Retorna B_lambda en unidades SI (W sr^-1 m^-3).
    l_m puede ser numpy array.
    """
    l = np.array(l_m, dtype=float)
    # evitar división por 0
    with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
        x = (h * c) / (l * kB * T)
        denom = np.expm1(x)  # exp(x) - 1 de forma estable
        B = (2 * h * c**2) / (l**5 * denom)
    B = np.where(np.isfinite(B), B, 0.0)
    return B

# ----------------------- Muestreo y ruido ----------------------------------
def sample_spectrum(T, l_min_um, l_max_um, N, pattern='uniform'):
    """
    Genera muestras en micrómetros y metros y la radiancia teórica (sin ruido).
    pattern: 'uniform', 'log', 'peak_concentrated'
    """
    if pattern == 'uniform':
        l_um = np.linspace(l_min_um, l_max_um, N)
    elif pattern == 'log':
        l_um = np.exp(np.linspace(np.log(l_min_um), np.log(l_max_um), N))
    elif pattern == 'peak_concentrated':
        l_max_theo_um = WIEN_B / T * 1e6
        n_peak = max(2, N // 2)
        n_rest = N - n_peak
        peak_range = (0.8 * l_max_theo_um, 1.2 * l_max_theo_um)
        peak_points = np.linspace(peak_range[0], peak_range[1], n_peak)
        rest_points = np.linspace(l_min_um, l_max_um, n_rest)
        l_um = np.sort(np.concatenate([peak_points, rest_points]))
    else:
        raise ValueError("pattern debe ser 'uniform', 'log' o 'peak_concentrated'")
    l_m = l_um * 1e-6
    I_true = planck_lambda(l_m, T)
    return l_um, l_m, I_true

def add_noise(I_true, noise_type='gaussian', rel_sigma=0.05, seed=None):
    rng = np.random.default_rng(seed)
    if noise_type == 'gaussian':
        sigma = rel_sigma * np.max(I_true)
        noise = rng.normal(0.0, sigma, size=I_true.shape)
        I_obs = I_true + noise
        I_obs = np.where(I_obs > 0, I_obs, 0.0)
        return I_obs
    elif noise_type == 'poisson':
        scale = 1e6 / np.max(I_true)
        counts = rng.poisson(I_true * scale)
        return counts / scale
    else:
        raise ValueError("noise_type debe ser 'gaussian' o 'poisson'")

# ----------------------- Reconstrucción ------------------------------------
def reconstruct_cubic(lambda_m, I_obs):
    return CubicSpline(lambda_m, I_obs, bc_type='natural')

def reconstruct_smoothing(lambda_m, I_obs, s=None, max_attempts=10):
    """
    Construye UnivariateSpline con parámetro s robusto.
    Si el algoritmo falla por 's too small', aumentamos s iterativamente.
    """
    # Estimación inicial si no se dio s
    N = len(lambda_m)
    if s is None:
        sigma_est = (np.max(I_obs) - np.min(I_obs)) * 0.01  # estimación conservadora
        s = (sigma_est**2) * N

    attempt = 0
    current_s = float(s)
    while attempt < max_attempts:
        try:
            # Suprimir warnings temporales (el wrapper de scipy puede emitir advertencias)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                sps = UnivariateSpline(lambda_m, I_obs, s=current_s)
            # Si construyó sin excepción, devolvemos
            return sps, current_s
        except Exception as e:
            # En caso de fallo, aumentar s y reintentar
            attempt += 1
            current_s *= 10.0
    # Si aún falla, levantar excepción con mensaje claro
    raise RuntimeError(f"Imposible construir smoothing spline tras {max_attempts} intentos. Último s={current_s}")

def reconstruct_planck_fit(lambda_m, I_obs, scale_guess=None, T_guess=None, bounds=None):
    """
    Ajuste I_obs ≈ A * B_lambda(lambda, T).
    bounds: tuple (lower_array, upper_array) o None
    """
    if T_guess is None:
        idx_max = np.nanargmax(I_obs)
        l_guess = lambda_m[idx_max]
        T_guess = WIEN_B / (l_guess if l_guess > 0 else 1e-6)

    if scale_guess is None:
        scale_guess = np.max(I_obs) / (np.max(planck_lambda(lambda_m, T_guess)) + 1e-30)

    def model(l, T, A):
        return A * planck_lambda(l, T)

    if bounds is None:
        # límites razonables: T entre 1K y 1e6 K, A entre 0 y inf
        lower = [1.0, 0.0]
        upper = [1e7, np.inf]
        bounds = (lower, upper)

    popt, pcov = curve_fit(model, lambda_m, I_obs, p0=[T_guess, scale_guess], bounds=bounds, maxfev=100000)
    T_fit, A_fit = popt
    def I_fit(l):
        return A_fit * planck_lambda(l, T_fit)
    return I_fit, (T_fit, A_fit), pcov

# ----------------------- Máximo y área -------------------------------------
def find_lambda_max_from_spline(spline, lam_min, lam_max, num_intervals=1200):
    """
    Encuentra lambda_max resolviendo s'(lambda)=0 con muestreo y brentq.
    Usa spline.derivative(k) (posicional).
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
            except ValueError:
                pass
    # Filtrar raíces y elegir la que tenga s'' < 0 y mayor valor de s
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
    # fallback: máximo de la malla si no se encontraron raíces
    lam_fine = np.linspace(lam_min, lam_max, max(1000, num_intervals))
    vals = spline(lam_fine)
    idx = np.nanargmax(vals)
    return lam_fine[idx]

def integrate_trapz(fun_or_spline, a, b, num=1200):
    x = np.linspace(a, b, num)
    if callable(fun_or_spline):
        y = fun_or_spline(x)
    else:
        y = fun_or_spline(x)
    mask = np.isfinite(y)
    if np.count_nonzero(mask) < 2:
        return 0.0
    return np.trapz(y[mask], x[mask])

# ----------------------- Bootstrap -----------------------------------------
def bootstrap_estimates(lambda_m, I_obs, method='smoothing', s=None, nboot=1500, seed=0):
    rng = np.random.default_rng(seed)
    lam = np.array(lambda_m)
    I = np.array(I_obs)
    n = len(lam)
    lam_min, lam_max = lam.min(), lam.max()

    est_lambda_max = []
    est_T_wien = []
    est_area = []

    for _ in range(nboot):
        idx = rng.integers(0, n, n)
        lam_b = lam[idx]
        I_b = I[idx]
        order = np.argsort(lam_b)
        lam_b = lam_b[order]
        I_b = I_b[order]
        try:
            if method == 'cubic':
                sc = reconstruct_cubic(lam_b, I_b)
                lammax = find_lambda_max_from_spline(sc, lam_min, lam_max)
                area = integrate_trapz(sc, lam_min, lam_max)
            elif method == 'smoothing':
                sps, _ = reconstruct_smoothing(lam_b, I_b, s=s)
                lammax = find_lambda_max_from_spline(sps, lam_min, lam_max)
                area = integrate_trapz(sps, lam_min, lam_max)
            elif method == 'planck':
                Ifit, params, _ = reconstruct_planck_fit(lam_b, I_b)
                Tfit, _ = params
                lammax = WIEN_B / Tfit
                area = integrate_trapz(Ifit, lam_min, lam_max)
            else:
                raise ValueError("method must be 'cubic','smoothing' or 'planck'")
            est_lambda_max.append(lammax)
            est_T_wien.append(WIEN_B / lammax)
            est_area.append(area)
        except Exception:
            # ignorar muestra bootstrap si falla
            continue

    return (np.array(est_lambda_max), np.array(est_T_wien), np.array(est_area))

def bootstrap_stats(b):
    return {
        'mean': np.nanmean(b),
        'std': np.nanstd(b, ddof=1),
        'pinf': np.nanpercentile(b, 2.5),
        'psup': np.nanpercentile(b, 97.5)
    }

# ----------------------- Experimento completo ------------------------------
def run_experiment(T, lmin_um, lmax_um, pattern, N,
                   noise_alpha, noise_type='gaussian',
                   s_smoothing=None, nboot=100, save_prefix="exp"):
    """
    Ejecuta un experimento (genera datos, reconstruye, estima, guarda).
    """
    # 1) datos
    l_um, l_m, I_true = sample_spectrum(T, lmin_um, lmax_um, N, pattern=pattern)
    I_obs = add_noise(I_true, noise_type=noise_type, rel_sigma=noise_alpha, seed=42)
    lam = l_m
    lam_min, lam_max = lam.min(), lam.max()

    # 2) reconstrucciones
    # cubic
    try:
        spline_cubic = reconstruct_cubic(lam, I_obs)
    except Exception as e:
        print("Error en CubicSpline:", e)
        spline_cubic = None

    # smoothing (robusta)
    try:
        sps_obj = None
        s_used = None
        sps_obj, s_used = reconstruct_smoothing(lam, I_obs, s=s_smoothing)
    except RuntimeError as e:
        # si falla, intentar con valor mayor de s incrementalmente (ya maneja internamente)
        print("Warning: smoothing spline inicial falló:", e)
        # reintentar con s mayor (approx)
        try:
            sps_obj, s_used = reconstruct_smoothing(lam, I_obs, s=((np.max(I_obs) - np.min(I_obs))**2) * N * 10)
        except Exception as e2:
            print("Smoothing spline falló nuevamente:", e2)
            sps_obj = None
            s_used = None
    except Exception as e:
        print("Error inesperado en smoothing spline:", e)
        sps_obj = None
        s_used = None

    # planck fit
    try:
        I_planck_fit_func, params_planck, pcov = reconstruct_planck_fit(lam, I_obs)
        T_fit_planck = params_planck[0]
    except Exception as e:
        print("Ajuste Planck falló:", e)
        I_planck_fit_func = None
        params_planck = (np.nan, np.nan)
        T_fit_planck = np.nan

    # 3) estimaciones
    lmax_true = WIEN_B / T
    true_area_partial = integrate_trapz(lambda x: planck_lambda(x, T), lam_min, lam_max)

    # cubic
    if spline_cubic is not None:
        lammax_cubic = find_lambda_max_from_spline(spline_cubic, lam_min, lam_max)
        T_cubic_wien = WIEN_B / lammax_cubic
        area_cubic = integrate_trapz(spline_cubic, lam_min, lam_max)
    else:
        lammax_cubic = np.nan; T_cubic_wien = np.nan; area_cubic = np.nan

    # smoothing
    if sps_obj is not None:
        lammax_smooth = find_lambda_max_from_spline(sps_obj, lam_min, lam_max)
        T_smooth_wien = WIEN_B / lammax_smooth
        area_smooth = integrate_trapz(sps_obj, lam_min, lam_max)
    else:
        lammax_smooth = np.nan; T_smooth_wien = np.nan; area_smooth = np.nan

    # planck
    if I_planck_fit_func is not None and isfinite(T_fit_planck):
        lammax_planck = WIEN_B / T_fit_planck
        area_planck = integrate_trapz(I_planck_fit_func, lam_min, lam_max)
    else:
        lammax_planck = np.nan; area_planck = np.nan

    # 4) bootstrap (incertidumbre)
    b_cubic = bootstrap_estimates(lam, I_obs, method='cubic', s=None, nboot=nboot, seed=123)
    b_smooth = bootstrap_estimates(lam, I_obs, method='smoothing', s=s_used, nboot=nboot, seed=123)
    b_planck = bootstrap_estimates(lam, I_obs, method='planck', s=None, nboot=nboot, seed=123)

    bs_cubic = bootstrap_stats(b_cubic)
    bs_smooth = bootstrap_stats(b_smooth)
    bs_planck = bootstrap_stats(b_planck)

    # 5) figura
    lam_fine = np.linspace(lam_min, lam_max, 800)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.scatter(l_um, I_obs, label='Datos (ruidosos)', s=18)
    ax.plot(l_um, I_true, label='Teórico (Planck)', linewidth=1)
    if spline_cubic is not None:
        ax.plot(lam_fine * 1e6, spline_cubic(lam_fine), label='Spline cúbico (interpol.)', linewidth=1)
    if sps_obj is not None:
        ax.plot(lam_fine * 1e6, sps_obj(lam_fine), label=f'Spline suavizado (s≈{s_used:.3e})', linewidth=1)
    if I_planck_fit_func is not None:
        ax.plot(lam_fine * 1e6, I_planck_fit_func(lam_fine), label=f'Ajuste Planck (T≈{T_fit_planck:.1f} K)', linewidth=1)
    ax.axvline(lmax_true * 1e6, linestyle='--', label='λ_max verdadero\n λ_max cubic\nλ_max smooth\nλ_max Planck-fit', alpha=0.6)
    ax.set_xlabel('λ (µm)')
    ax.set_ylabel('Radiancia B_λ (SI)')
    ax.set_title(f'Reconstrucción: T={T} K, N={N}, ruido={noise_alpha*100:.1f}% , patrón={pattern}')
    ax.legend(fontsize='small', loc='upper right')
    fig.tight_layout()
    fig_name = f"{save_prefix}_recon.png"
    fig_path = save_fig_safe(fig, fig_name)

    # 6) resumen CSV con experimentos en columnas - TODOS EN UN ARCHIVO
    summary = {
        'experiment': save_prefix,
        'T_true': T,
        'N': N,
        'pattern': pattern,
        'noise_alpha': noise_alpha,
        'lam_max_true_um': lmax_true * 1e6,
        'lam_max_cubic_um': lammax_cubic * 1e6 if isfinite(lammax_cubic) else np.nan,
        'lam_max_smooth_um': lammax_smooth * 1e6 if isfinite(lammax_smooth) else np.nan,
        'lam_max_planck_um': lammax_planck * 1e6 if isfinite(lammax_planck) else np.nan,
        'T_cubic_wien': T_cubic_wien,
        'T_smooth_wien': T_smooth_wien,
        'T_planck_fit': T_fit_planck,
        'area_true_partial': true_area_partial,
        'area_cubic': area_cubic,
        'area_smooth': area_smooth,
        'area_planck': area_planck,

        # --- BOOTSTRAP: CUBIC ---
        'bs_cubic_mean': bs_cubic['mean'],
        'bs_cubic_std': bs_cubic['std'],
        'bs_cubic_pinf': bs_cubic['pinf'],
        'bs_cubic_psup': bs_cubic['psup'],

        # --- BOOTSTRAP: SMOOTH ---
        'bs_smooth_mean': bs_smooth['mean'],
        'bs_smooth_std': bs_smooth['std'],
        'bs_smooth_pinf': bs_smooth['pinf'],
        'bs_smooth_psup': bs_smooth['psup'],

        # --- BOOTSTRAP: PLANCK ---
        'bs_planck_mean': bs_planck['mean'],
        'bs_planck_std': bs_planck['std'],
        'bs_planck_pinf': bs_planck['pinf'],
        'bs_planck_psup': bs_planck['psup'],
    }

    # CAMBIO CRÍTICO: Usar un nombre de archivo común para todos los experimentos
    csv_name = "all_experiments_summary.csv"  # <-- ARCHIVO ÚNICO
    csv_path = os.path.join(RESULTS_DIR, csv_name)

    # Crear DataFrame vertical: parámetros en primera columna, valores en columna del experimento
    df_new = pd.DataFrame({
        'Parameter': list(summary.keys()),
        save_prefix: list(summary.values())
    })
    df_new.set_index('Parameter', inplace=True)

    # Si el archivo existe, agregar como nueva columna
    if os.path.exists(csv_path):
        df_old = pd.read_csv(csv_path, index_col=0)
        df_combined = pd.concat([df_old, df_new], axis=1)
        df_combined.to_csv(csv_path, quoting=1)  # quoting=1 protege las comas
    else:
        df_new.to_csv(csv_path, quoting=1)

    # guardar bootstrap resultados
    np.savez(os.path.join(RESULTS_DIR, f"{save_prefix}_bootstrap.npz"),
             cubic=b_cubic, smooth=b_smooth, planck=b_planck)

    return {
        'summary': summary,
        'fig_path': fig_path,
        'csv_path': csv_path,
        'bootstrap': {'cubic': b_cubic, 'smooth': b_smooth, 'planck': b_planck}
    }

# ----------------------- Ejecución de experimentos mínimos -----------------
if __name__ == "__main__":
    # Experimento A (principal)
    print("Ejecutando Experimento A (T=5800 K, N=20, ruido 5%) ...")
    resA = run_experiment(T=5800, lmin_um=0.1, lmax_um=3.0, pattern='uniform',
                          N=20, noise_alpha=0.05, noise_type='gaussian',
                          s_smoothing=None, nboot=80, save_prefix="expA")
    print("Guardado:", resA['fig_path'], resA['csv_path'])

    # Experimento B (sensibilidad al muestreo)
    for N in [10, 30, 100]:
        print("Ejecutando Experimento B (N=", N, ") ...")
        resB = run_experiment(T=5800, lmin_um=0.1, lmax_um=3.0, pattern='uniform',
                              N=N, noise_alpha=0.05, noise_type='gaussian',
                              s_smoothing=None, nboot=40, save_prefix=f"expB_N{N}")
        print("Guardado:", resB['fig_path'], resB['csv_path'])

    # Experimento C (sensibilidad al ruido)
    for alpha in [0.01, 0.05, 0.10]:
        print("Ejecutando Experimento C (ruido=", alpha, ") ...")
        resC = run_experiment(T=3000, lmin_um=0.3, lmax_um=5.0, pattern='uniform',
                              N=30, noise_alpha=alpha, noise_type='gaussian',
                              s_smoothing=None, nboot=40, save_prefix=f"expC_noise{int(alpha*100)}")
        print("Guardado:", resC['fig_path'], resC['csv_path'])

    print("Todos los experimentos mínimos completados. Revisa la carpeta 'results' junto al script.")
