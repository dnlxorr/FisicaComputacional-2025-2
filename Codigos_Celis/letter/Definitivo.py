import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson
from numpy.polynomial.legendre import leggauss


# ============================================================
# 1. PARÁMETROS FÍSICOS
# ============================================================

G = 6.674e-11
M = 5.972e24
mu = G*M

a = 7e6
e = 0.2


# ============================================================
# 2. ORBITA REAL
# ============================================================

def orbit_true(n_points=2000):
    theta = np.linspace(0, 2*np.pi, n_points)
    r = a*(1 - e**2)/(1 - e*np.cos(theta))
    x = r*np.cos(theta)
    y = r*np.sin(theta)
    return theta, x, y


# ============================================================
# 3. INTERPOLACIÓN DE LAGRANGE (GENÉRICA 1D)
# ============================================================

def lagrange_interpolation(x_nodes, y_nodes, x_eval):
    L = np.zeros_like(x_eval)
    n = len(x_nodes)

    for i in range(n):
        li = np.ones_like(x_eval)
        for j in range(n):
            if i != j:
                li *= (x_eval - x_nodes[j])/(x_nodes[i] - x_nodes[j])
        L += y_nodes[i]*li
    return L


# ============================================================
# 4. EXPERIMENTO 1 — VARIAR NODOS
# ============================================================

def experiment_n_nodes(n_list=[5, 8, 12, 16, 20, 25]):
    theta, x_true, y_true = orbit_true()
    theta_eval = np.linspace(0, 2*np.pi, 4000)

    results = []

    for N in n_list:

        # Selección de N nodos
        idx = np.linspace(0, len(theta)-1, N, dtype=int)
        th_nodes = theta[idx]
        x_nodes = x_true[idx]
        y_nodes = y_true[idx]

        # Interpolación
        x_interp = lagrange_interpolation(th_nodes, x_nodes, theta_eval)
        y_interp = lagrange_interpolation(th_nodes, y_nodes, theta_eval)

        # Órbita real evaluada
        x_true_eval = np.interp(theta_eval, theta, x_true)
        y_true_eval = np.interp(theta_eval, theta, y_true)

        # Error euclidiano
        error = np.sqrt((x_interp - x_true_eval)**2 +
                         (y_interp - y_true_eval)**2)
        error_max = np.max(error)
        results.append((N, error_max))

        # ====================================================
        # FIGURA COMBINADA
        # ====================================================
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        # --- Panel 1: Órbita ---
        axes[0].plot(x_true_eval, y_true_eval, label="Órbita real", linewidth=2)
        axes[0].plot(x_interp, y_interp, "--", label=f"Interpolada (N={N})")
        axes[0].set_aspect("equal")
        axes[0].set_xlabel("x (m)")
        axes[0].set_ylabel("y (m)")
        axes[0].set_title(f"Órbita real vs interpolada (N={N})")
        axes[0].grid(True)
        axes[0].legend()

        # --- Panel 2: Error ---
        axes[1].plot(theta_eval, error)
        axes[1].set_xlabel(r"$\theta$")
        axes[1].set_ylabel("Error (m)")
        axes[1].set_title(f"Error euclidiano (máx = {error_max:.2e})")
        axes[1].grid(True)

        plt.tight_layout()
        plt.show()

    return results



# ============================================================
# 5. EXPERIMENTO 2 — VARIAR EXCENTRICIDAD
# ============================================================

def experiment_eccentricity(e_list=[0.0, 0.2, 0.4, 0.6, 0.8]):

    N = 15
    theta_eval = np.linspace(0, 2*np.pi, 4000)
    errors = []

    # Crear figura con subplots
    fig, axes = plt.subplots(1, len(e_list), figsize=(4*len(e_list), 4))

    for k, ei in enumerate(e_list):
        global e
        e = ei

        # Órbita real
        theta, x_true, y_true = orbit_true()
        x_true_eval = np.interp(theta_eval, theta, x_true)
        y_true_eval = np.interp(theta_eval, theta, y_true)

        # Nodos
        idx = np.linspace(0, len(theta)-1, N, dtype=int)
        th_nodes = theta[idx]
        x_nodes = x_true[idx]
        y_nodes = y_true[idx]

        # Interpolación
        x_interp = lagrange_interpolation(th_nodes, x_nodes, theta_eval)
        y_interp = lagrange_interpolation(th_nodes, y_nodes, theta_eval)

        # Error máximo
        error = np.max(np.sqrt((x_interp - x_true_eval)**2 +
                               (y_interp - y_true_eval)**2))
        errors.append((ei, error))

        # ---- Gráfica ----
        ax = axes[k]
        ax.plot(x_true_eval, y_true_eval, label="Real", linewidth=2)
        ax.plot(x_interp, y_interp, "--", label="Interpolada")
        ax.set_aspect("equal")
        ax.set_title(rf"$e={ei}$" + f"\nError máx = {error:.2e}")
        ax.grid(True)

        if k == 0:
            ax.set_ylabel("y (m)")
        ax.set_xlabel("x (m)")

    axes[0].legend(loc="upper right")
    plt.suptitle("Interpolación orbital para distintas excentricidades", fontsize=14)
    plt.tight_layout()
    plt.show()

    return errors



# ============================================================
# 6. EXPERIMENTO 3 — INTEGRACIÓN NUMÉRICA (REAL + INTERPOLADA)
# ============================================================

def experiment_integration():

    # ----- ÓRBITA REAL -----
    theta, x_true, y_true = orbit_true()
    r_true = np.sqrt(x_true**2 + y_true**2)

    # Integral real (Simpson)
    I_simpson_real = simpson(r_true, theta)

    # Integral real (Gauss-Legendre)
    xg, wg = leggauss(8)
    mid = 0.5*(theta[-1] + theta[0])
    s = 0.5*(theta[-1] - theta[0])
    I_gauss_real = s * np.sum(wg * np.interp(mid + s*xg, theta, r_true))

    # ----- INTERPOLACIÓN -----
    N = 20
    idx = np.linspace(0, len(theta)-1, N, dtype=int)
    th_nodes = theta[idx]

    r_nodes = r_true[idx]
    r_interp_eval = lagrange_interpolation(th_nodes, r_nodes, theta)

    # Integral interpolada (Simpson)
    I_simpson_interp = simpson(r_interp_eval, theta)

    # Integral interpolada (Gauss-Legendre)
    I_gauss_interp = s * np.sum(wg * np.interp(mid + s*xg, theta, r_interp_eval))

    # ----- RESULTADOS -----
    print("\n============================")
    print("   INTEGRACIÓN COMPARADA    ")
    print("============================\n")

    print(">>> Integral real:")
    print("Simpson real        =", I_simpson_real)
    print("Gauss-Legendre real =", I_gauss_real)

    print("\n>>> Integral interpolada:")
    print("Simpson interpolada        =", I_simpson_interp)
    print("Gauss-Legendre interpolada =", I_gauss_interp)

    print("\n>>> Diferencias:")
    print("Simpson: |real - interp| =", abs(I_simpson_real - I_simpson_interp))
    print("Gauss:   |real - interp| =", abs(I_gauss_real - I_gauss_interp))


def experiment_fixed_e_vary_nodes(
    e_fixed=0.8,
    n_list=[6, 10, 14, 18, 22, 26]
):
    global e
    e = e_fixed

    theta, x_true, y_true = orbit_true()
    theta_eval = np.linspace(0, 2*np.pi, 4000)

    x_true_eval = np.interp(theta_eval, theta, x_true)
    y_true_eval = np.interp(theta_eval, theta, y_true)

    # Figura con subgráficas
    fig, axes = plt.subplots(1, len(n_list), figsize=(4*len(n_list), 4))

    for k, N in enumerate(n_list):

        # Selección de nodos
        idx = np.linspace(0, len(theta)-1, N, dtype=int)
        th_nodes = theta[idx]
        x_nodes = x_true[idx]
        y_nodes = y_true[idx]

        # Interpolación
        x_interp = lagrange_interpolation(th_nodes, x_nodes, theta_eval)
        y_interp = lagrange_interpolation(th_nodes, y_nodes, theta_eval)

        # Error máximo
        error = np.max(np.sqrt((x_interp - x_true_eval)**2 +
                               (y_interp - y_true_eval)**2))

        # Gráfica
        ax = axes[k]
        ax.plot(x_true_eval, y_true_eval, label="Órbita real", linewidth=2)
        ax.plot(x_interp, y_interp, "--", label="Interpolada")
        ax.set_aspect("equal")
        ax.set_title(f"N={N}\nError máx = {error:.2e}")
        ax.grid(True)

        if k == 0:
            ax.set_ylabel("y (m)")
        ax.set_xlabel("x (m)")

    axes[0].legend(loc="upper right")
    plt.suptitle(
        rf"Ajuste de la interpolación para excentricidad fija $e={e_fixed}$",
        fontsize=14
    )
    plt.tight_layout()
    plt.show()

# ============================================================
# 7. MAIN
# ============================================================

if __name__ == "__main__":

    print("\n=== EXPERIMENTO 1: VARIAR NODOS ===")
    res1 = experiment_n_nodes()
    print(res1)

    print("\n=== EXPERIMENTO 2: VARIAR EXCENTRICIDAD ===")
    res2 = experiment_eccentricity()
    print(res2)

    print("\n=== EXPERIMENTO 3: INTEGRACIÓN ===")
    experiment_integration()

    print("\n=== EXPERIMENTO 4: EXCENTRICIDAD FIJA, NODOS VARIABLES ===")
    experiment_fixed_e_vary_nodes()

