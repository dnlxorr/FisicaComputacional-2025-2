# Python code to reconstruct blackbody spectrum using splines, compute integrals and derivatives,
# and produce diagnostic plots and error metrics.
# This will run in the notebook environment and display plots and numeric results.
# Requires: numpy, scipy, matplotlib, pandas (these are usually available).
# If any import fails here, tell me and I'll adapt the code.

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import h, c, k, sigma
from scipy.interpolate import CubicSpline, UnivariateSpline
from scipy.optimize import brentq
import pandas as pd

plt.rcParams['figure.dpi'] = 100


def planck_lambda(lam, T):
    """Planck spectral radiance per unit wavelength B(λ,T) [W m^-2 sr^-1 m^-1].
    lam and T in SI (m, K)."""
    # Avoid overflow for very small lambda by using safe exponent
    a = 2 * h * c ** 2 / (lam ** 5)

    b = h * c / (lam * k * T)
    # handle large b safely
    with np.errstate(over='ignore', divide='ignore'):
        return a / (np.expm1(b))


def analytic_integral_B_over_lambda(T):
    """Integral from 0 to inf of B(λ,T) dλ = (sigma T^4) / pi"""
    return sigma * T ** 4 / np.pi


def find_spline_peak(cs, lam_grid):
    """Find lambda that maximizes the cubic spline cs on interval defined by lam_grid.
    We search for stationary points per interval by solving derivative = 0 for each cubic segment.
    cs is scipy.interpolate.CubicSpline (piecewise cubic with .c attribute)."""
    # Use derivative spline
    dcs = cs.derivative()
    peaks = []
    # Look for roots of derivative in each interval by scanning sign changes and solving
    for i in range(len(lam_grid) - 1):
        a = lam_grid[i]
        b = lam_grid[i + 1]
        fa = dcs(a)
        fb = dcs(b)
        if np.isnan(fa) or np.isnan(fb):
            continue
        # If sign change, find root
        if fa == 0:
            peaks.append(a)
        if fa * fb < 0:
            try:
                root = brentq(lambda x: dcs(x), a, b)
                peaks.append(root)
            except Exception:
                pass
    # Also include endpoints
    peaks.extend([lam_grid[0], lam_grid[-1]])
    # Evaluate spline at candidates and pick max
    if len(peaks) == 0:
        return None, None
    peaks = np.array(peaks)
    vals = cs(peaks)
    idx = np.nanargmax(vals)
    return peaks[idx], vals[idx]


def reconstruct_and_evaluate(T=5800.0, N_samples=20, lam_min=1e-9, lam_max=3e-6,
                             sample_distribution='linear', noise_std_frac=0.0,
                             smoothing_spline_s=None):
    """Generate spectrum, sample, build splines (interpolating and smoothing), evaluate integrals and peaks."""
    # high-resolution reference grid for plotting and "truth" integration
    lam_ref = np.linspace(lam_min, lam_max, 5000)
    B_ref = planck_lambda(lam_ref, T)

    # choose sample points
    if sample_distribution == 'linear':
        lam_samp = np.linspace(lam_min, lam_max, N_samples)
    elif sample_distribution == 'log':
        lam_samp = np.exp(np.linspace(np.log(lam_min), np.log(lam_max), N_samples))
    elif sample_distribution == 'random':
        lam_samp = np.sort(np.random.uniform(lam_min, lam_max, N_samples))
    else:
        raise ValueError("Unknown distribution")

    B_samp = planck_lambda(lam_samp, T)
    # add gaussian noise relative to local spectral value
    if noise_std_frac > 0.0:
        noise = np.random.normal(scale=noise_std_frac * np.maximum(B_samp, 1e-30))
        B_noisy = B_samp + noise
    else:
        B_noisy = B_samp.copy()

    # Interpolating cubic spline (natural boundary conditions by default in CubicSpline with bc_type='natural')
    cs_interp = CubicSpline(lam_samp, B_noisy, bc_type='natural', extrapolate=False)
    B_cs_ref = cs_interp(lam_ref)

    # Smoothing spline (UnivariateSpline) - it fits and has smoothing factor s (if None it interpolates)
    if smoothing_spline_s is None:
        # default small smoothing (interpolating)
        smoothing_spline_s = 0.0
    uspline = UnivariateSpline(lam_samp, B_noisy, s=smoothing_spline_s, k=3)
    B_us_ref = uspline(lam_ref)

    # Integrals: exact analytic for B_ref (over domain lam_min..lam_max) using numerical integrate since full analytic from 0..inf known
    # We'll compute:
    #  - reference integral of exact Planck over [lam_min, lam_max] by numerical integration (trapz on fine grid)
    #  - spline integral using antiderivative for CubicSpline and UnivariateSpline.integral for smoothing spline
    integral_ref = np.trapz(B_ref, lam_ref)  # integral over limited domain
    total_true = analytic_integral_B_over_lambda(T)  # integral 0..inf
    # CubicSpline antiderivative: use .antiderivative()
    cs_ant = cs_interp.antiderivative()
    try:
        integral_cs = cs_ant(lam_max) - cs_ant(lam_min)
    except Exception:
        integral_cs = np.nan
    # UnivariateSpline integral method
    try:
        integral_us = uspline.integral(lam_min, lam_max)
    except Exception:
        integral_us = np.nan

    # Peak locations
    # For the true Planck on lam_ref grid, approximate peak
    idx_peak_ref = np.nanargmax(B_ref)
    lam_peak_ref = lam_ref[idx_peak_ref]
    B_peak_ref = B_ref[idx_peak_ref]
    # For cubic spline use find_spline_peak
    lam_peak_cs, B_peak_cs = find_spline_peak(cs_interp, lam_samp)
    # For smoothing spline we can find maximum by scanning lam_ref
    idx_peak_us = np.nanargmax(B_us_ref)
    lam_peak_us = lam_ref[idx_peak_us]
    B_peak_us = B_us_ref[idx_peak_us]

    # RMSE on lam_ref grid (compare reconstructed to true)
    rmse_cs = np.sqrt(np.nanmean((B_cs_ref - B_ref) ** 2))
    rmse_us = np.sqrt(np.nanmean((B_us_ref - B_ref) ** 2))

    # Relative integral error vs analytic total_true (0..inf)
    rel_err_cs_vs_total = (integral_cs - total_true) / total_true if np.isfinite(integral_cs) else np.nan
    rel_err_us_vs_total = (integral_us - total_true) / total_true if np.isfinite(integral_us) else np.nan
    rel_err_cs_vs_refdomain = (integral_cs - integral_ref) / (integral_ref if integral_ref != 0 else np.nan)
    rel_err_us_vs_refdomain = (integral_us - integral_ref) / (integral_ref if integral_ref != 0 else np.nan)

    # Pack results
    results = {
        'lam_ref': lam_ref, 'B_ref': B_ref,
        'lam_samp': lam_samp, 'B_samp': B_samp, 'B_noisy': B_noisy,
        'cs': cs_interp, 'us': uspline,
        'B_cs_ref': B_cs_ref, 'B_us_ref': B_us_ref,
        'integral_ref_domain': integral_ref, 'integral_cs': integral_cs, 'integral_us': integral_us,
        'total_true_0_inf': total_true,
        'lam_peak_ref': lam_peak_ref, 'B_peak_ref': B_peak_ref,
        'lam_peak_cs': lam_peak_cs, 'B_peak_cs': B_peak_cs,
        'lam_peak_us': lam_peak_us, 'B_peak_us': B_peak_us,
        'rmse_cs': rmse_cs, 'rmse_us': rmse_us,
        'rel_err_cs_vs_total': rel_err_cs_vs_total, 'rel_err_us_vs_total': rel_err_us_vs_total,
        'rel_err_cs_vs_refdomain': rel_err_cs_vs_refdomain, 'rel_err_us_vs_refdomain': rel_err_us_vs_refdomain,
    }
    return results


# Example run with default parameters (you can modify below and re-run)
T = 5800.0  # temperature in K
N_samples = 20
lam_min = 1e-9
lam_max = 3e-6
sample_distribution = 'linear'  # 'linear', 'log', or 'random'
noise_std_frac = 0.02  # relative noise (2%); set 0 to remove noise
smoothing_spline_s = 1e-26  # smoothing factor for UnivariateSpline; 0 -> interpolation

res = reconstruct_and_evaluate(T=T, N_samples=N_samples, lam_min=lam_min, lam_max=lam_max,
                               sample_distribution=sample_distribution,
                               noise_std_frac=noise_std_frac,
                               smoothing_spline_s=smoothing_spline_s)

# Print numeric summary
summary = {
    "Temperature (K)": T,
    "Samples (N)": N_samples,
    "Sampling dist.": sample_distribution,
    "Noise std frac": noise_std_frac,
    "Analytic integral 0..inf (B dλ)": res['total_true_0_inf'],
    "Integral reference domain (trapz)": res['integral_ref_domain'],
    "Integral CubicSpline (lam_min..lam_max)": res['integral_cs'],
    "Integral SmoothingSpline (lam_min..lam_max)": res['integral_us'],
    "RMSE CubicSpline (vs ref grid)": res['rmse_cs'],
    "RMSE SmoothingSpline (vs ref grid)": res['rmse_us'],
    "Peak lam (ref grid)": res['lam_peak_ref'],
    "Peak lam (CubicSpline)": res['lam_peak_cs'],
    "Peak lam (SmoothingSpline)": res['lam_peak_us'],
    "Rel err CubicSpline vs total_true": res['rel_err_cs_vs_total'],
    "Rel err SmoothingSpline vs total_true": res['rel_err_us_vs_total'],
}

df_summary = pd.DataFrame(list(summary.items()), columns=['Quantity', 'Value'])
import caas_jupyter_tools as cjt;

cjt.display_dataframe_to_user("Resumen de resultados", df_summary)

# Plot 1: Spectra (reference, sampled, splines)
plt.figure(figsize=(8, 4))
plt.plot(res['lam_ref'] * 1e9, res['B_ref'], label='Planck (ref)')  # convert to nm for x-axis
plt.scatter(res['lam_samp'] * 1e9, res['B_noisy'], label='Samples (noisy)', zorder=10)
plt.plot(res['lam_ref'] * 1e9, res['B_cs_ref'], label='Cubic Spline interp')
plt.plot(res['lam_ref'] * 1e9, res['B_us_ref'], label='Smoothing Spline')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Spectral radiance B(λ,T) [W m$^{-2}$ sr$^{-1}$ m$^{-1}$]')
plt.title(f'Blackbody spectrum reconstruction, T={T} K')
plt.legend()
plt.tight_layout()
plt.show()

# Plot 2: Zoom around peak
plt.figure(figsize=(8, 4))
mask = (res['lam_ref'] >= res['lam_peak_ref'] * 0.3) & (res['lam_ref'] <= res['lam_peak_ref'] * 3.0)
plt.plot(res['lam_ref'][mask] * 1e9, res['B_ref'][mask], label='Planck (ref)')
plt.scatter(res['lam_samp'] * 1e9, res['B_noisy'], label='Samples (noisy)', zorder=10)
plt.plot(res['lam_ref'][mask] * 1e9, res['B_cs_ref'][mask], label='Cubic Spline')
plt.plot(res['lam_ref'][mask] * 1e9, res['B_us_ref'][mask], label='Smoothing Spline')
plt.axvline(res['lam_peak_ref'] * 1e9, linestyle='--', label='Peak ref')
if res['lam_peak_cs'] is not None:
    plt.axvline(res['lam_peak_cs'] * 1e9, linestyle=':', label='Peak cubic spline')
plt.axvline(res['lam_peak_us'] * 1e9, linestyle='-.', label='Peak smoothing spline')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Spectral radiance')
plt.title('Zoom near spectral peak')
plt.legend()
plt.tight_layout()
plt.show()

# Table of key numeric metrics
metrics = {
    "metric": ["integral_ref_domain", "integral_cs", "integral_us", "total_true_0_inf", "rmse_cs", "rmse_us",
               "lam_peak_ref_nm", "lam_peak_cs_nm", "lam_peak_us_nm", "rel_err_cs_vs_total", "rel_err_us_vs_total"],
    "value": [res['integral_ref_domain'], res['integral_cs'], res['integral_us'], res['total_true_0_inf'],
              res['rmse_cs'], res['rmse_us'],
              res['lam_peak_ref'] * 1e9, (res['lam_peak_cs'] * 1e9 if res['lam_peak_cs'] is not None else np.nan),
              (res['lam_peak_us'] * 1e9 if res['lam_peak_us'] is not None else np.nan),
              res['rel_err_cs_vs_total'], res['rel_err_us_vs_total']]
}
df_metrics = pd.DataFrame(metrics)
cjt.display_dataframe_to_user("Métricas numéricas clave", df_metrics)

print("He terminado la ejecución. Si quieres, puedo:")
print("- ejecutar más experimentos variando N_samples, ruido o distribución;")
print("- exportar las figuras como archivos PNG/SVG;")
print("- generar un script .py o un cuaderno Jupyter con este código;")
print("- integrar búsqueda de la mejor smoothing factor por validación cruzada;")
print("- añadir B-splines explícitos o integración sobre ángulos para obtener irradiancia total.")
