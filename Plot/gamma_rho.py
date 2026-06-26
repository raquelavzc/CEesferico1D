import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


def fit_line(x, y):
    slope, intercept = np.polyfit(x, y, 1)
    yfit = intercept + slope * x
    ss_res = np.sum((y - yfit) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return slope, intercept, r2, yfit


def main():
    parser = argparse.ArgumentParser(
        description="Grafica ln(rho_c_max) vs ln(phi_c - phi0) y ajusta gamma subcritica."
    )
    parser.add_argument("phi_c", type=float, help="Valor critico phi_c.")
    parser.add_argument("archivo", nargs="?", default="valores_centrales.dat")
    parser.add_argument("--output", default="gamma_rho_valores_centrales.png")
    parser.add_argument("--tabla", default="gamma_rho_valores_centrales.dat")
    parser.add_argument("--min-delta", type=float, default=None)
    parser.add_argument("--max-delta", type=float, default=None)
    args = parser.parse_args()

    data = np.genfromtxt(args.archivo, names=True)
    phi = data["phi0_real"]
    rho = data["rho_c_max"]
    delta = args.phi_c - phi

    mask = (delta > 0.0) & (rho > 0.0)
    if args.min_delta is not None:
        mask &= delta >= args.min_delta
    if args.max_delta is not None:
        mask &= delta <= args.max_delta

    phi = phi[mask]
    rho = rho[mask]
    delta = delta[mask]
    if len(phi) < 2:
        raise ValueError("Se necesitan al menos dos puntos validos para ajustar.")

    x = np.log(delta)
    y = np.log(rho)
    slope, intercept, r2, yfit = fit_line(x, y)
    gamma = -0.5 * slope

    order = np.argsort(x)
    fig, ax = plt.subplots(figsize=(8, 5), dpi=160)
    ax.scatter(x, y, color="#0b5d7a", s=45, label="Datos")
    ax.plot(
        x[order],
        yfit[order],
        color="#c2410c",
        lw=2,
        label=rf"Ajuste: $\gamma={gamma:.6f}$, $R^2={r2:.5f}$",
    )

    for xi, yi, di in zip(x, y, delta):
        ax.annotate(f"{di:.0e}", (xi, yi), textcoords="offset points", xytext=(4, 4), fontsize=8)

    ax.set_xlabel(r"$\ln(\phi_* - \phi_0)$")
    ax.set_ylabel(r"$\ln(\rho^{max}_{central})$")
    ax.set_title("Gamma subcritica por densidad central maxima")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.output)

    with Path(args.tabla).open("w", encoding="utf-8") as f:
        f.write("# Ajuste: ln(rho_c_max) = C - 2*gamma*ln(phi_c - phi0)\n")
        f.write(f"# phi_c = {args.phi_c:.17g}\n")
        f.write(f"# pendiente = {slope:.16e}\n")
        f.write(f"# gamma = {gamma:.16e}\n")
        f.write(f"# C = {intercept:.16e}\n")
        f.write(f"# R2 = {r2:.16e}\n")
        f.write("# phi0 rho_c_max delta ln_delta ln_rho\n")
        for row in zip(phi, rho, delta, x, y):
            f.write(" ".join(f"{value:.16e}" for value in row) + "\n")

    print(f"puntos usados = {len(phi)}")
    print(f"pendiente = {slope:.12g}")
    print(f"gamma = {gamma:.12g}")
    print(f"R2 = {r2:.12g}")
    print(f"grafica = {args.output}")
    print(f"tabla = {args.tabla}")


if __name__ == "__main__":
    main()
