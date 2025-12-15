import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline, UnivariateSpline

# Constantes de Planck
h = 6.62607015e-34  # J·s
c = 2.99792458e8    # m/s
k_B = 1.380649e-23  # J/K

def planck_lambda(wavelength_um, T):
    """Función de Planck en función de longitud de onda"""
    wavelength_m = wavelength_um * 1e-6
    factor = (2 * h * c**2) / wavelength_m**5
    exponent = (h * c) / (wavelength_m * k_B * T)
    return factor / (np.exp(exponent) - 1)

# Parámetros del problema
T = 5000  # Temperatura en K
lam_min, lam_max = 0.2, 2.5  # Rango de longitud de onda en µm

# Generar datos sintéticos (reemplazar con tus datos reales)
lambda_obs = np.linspace(lam_min, lam_max, 20)
I_obs = planck_lambda(lambda_obs, T) + np.random.normal(0, 1e13, len(lambda_obs))

# Crear splines
spline_cubic = CubicSpline(lambda_obs, I_obs)
spline_smooth = UnivariateSpline(lambda_obs, I_obs, s=1e28)  # s controla el suavizado

# Parámetros de ajuste (estos deberían venir de scipy.optimize.curve_fit)
T_fit = 4950
scale_fit = 1.05

# Lambda fina para curvas continuas
lam_fine = np.linspace(lam_min, lam_max, 2000)

# Curva teórica
I_planck_fine = planck_lambda(lam_fine, T)

# Reconstrucciones
I_cubic = spline_cubic(lam_fine)
I_smooth = spline_smooth(lam_fine)
I_planck_fit = scale_fit * planck_lambda(lam_fine, T_fit)

# Gráfico
plt.figure(figsize=(7, 4))
plt.plot(lam_fine, I_planck_fine, label='Planck teórico', linewidth=2)
plt.scatter(lambda_obs, I_obs, color='black', s=25, label='Datos discretos', zorder=5)
plt.plot(lam_fine, I_cubic, '--', label='Spline cúbico', alpha=0.8)
plt.plot(lam_fine, I_smooth, ':', label='Spline suavizado', linewidth=2)
plt.plot(lam_fine, I_planck_fit, '-.', label='Ajuste a Planck', alpha=0.8)

plt.xlabel('Longitud de onda [µm]')
plt.ylabel('Intensidad espectral (u.a.)')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('comparacion_espectros_A.png', dpi=300, bbox_inches='tight')
plt.show()