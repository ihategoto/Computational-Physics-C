#!/usr/bin/env python3
"""
Generate reference plots for all physlib examples.

Requirements: numpy, matplotlib
Run from anywhere — the script locates the project root automatically.
The examples must be built first (`make examples` from the project root).

Output: one subdirectory per example, each containing PNG files.
"""

import os, sys, subprocess
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Locate project root (one level above this script) ─────────────────────
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)

# ── Global style ──────────────────────────────────────────────────────────
plt.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 150,
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 12,
    "legend.fontsize": 10,
    "lines.linewidth": 1.4,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linestyle": "--",
})

TAB10 = plt.cm.tab10(np.linspace(0, 0.9, 10))


def run(binary, *args):
    """Run ./binary with the given arguments; return parsed stdout as (N, M) array."""
    cmd = ["./" + binary] + [str(a) for a in args]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"  [ERROR] {' '.join(cmd)}\n  {r.stderr.strip()}", file=sys.stderr)
        return np.empty((0, 0))
    rows = []
    for line in r.stdout.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        try:
            rows.append(list(map(float, line.split())))
        except ValueError:
            pass
    return np.array(rows) if rows else np.empty((0, 0))


def save(fig, path):
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"    {os.path.relpath(path)}")


# ══════════════════════════════════════════════════════════════════════════ #
# 1. Harmonic oscillator
# ══════════════════════════════════════════════════════════════════════════ #

def plot_harmonic_oscillator():
    out = "plots/harmonic_oscillator"
    os.makedirs(out, exist_ok=True)

    METHODS = ["euler", "euler_cromer", "leapfrog", "midpoint",
               "verlet_v", "pc", "rk2", "rk4"]
    LABELS  = {
        "euler":        "Euler",
        "euler_cromer": "Euler-Cromer",
        "leapfrog":     "Leapfrog",
        "midpoint":     "Midpoint",
        "verlet_v":     "Velocity Verlet",
        "pc":           "Pred.-Corr.",
        "rk2":          "RK2",
        "rk4":          "RK4",
    }
    SYMPLECTIC = {"euler_cromer", "leapfrog", "midpoint", "verlet_v"}
    COLORS = {m: TAB10[i] for i, m in enumerate(METHODS)}

    # Collect trajectories: t_max=50, dt=0.05, x0=1, v0=0, k=m=1
    data = {}
    for m in METHODS:
        d = run("harmonic_oscillator", 50, 0.05, 1.0, 0.0, 1.0, 1.0, m)
        if d.ndim == 2 and d.shape[1] >= 4:
            data[m] = d  # cols: t  x  v  (E-E0)/E0

    # ── Energy error (log scale) ──────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 5))
    for m in METHODS:
        if m not in data:
            continue
        t   = data[m][:, 0]
        err = np.abs(data[m][:, 3])
        err = np.where(err == 0, 1e-17, err)
        ls  = "-" if m in SYMPLECTIC else "--"
        ax.semilogy(t, err, ls, color=COLORS[m], label=LABELS[m])
    ax.set_xlabel("t")
    ax.set_ylabel(r"$|\Delta E / E_0|$")
    ax.set_title("Energy conservation — harmonic oscillator  (k = m = 1,  dt = 0.05)")
    ax.legend(ncol=2, framealpha=0.9)
    save(fig, f"{out}/energy_error.png")

    # ── Phase portraits (2 × 4 grid) ─────────────────────────────────────
    theta = np.linspace(0, 2 * np.pi, 400)
    fig, axes = plt.subplots(2, 4, figsize=(14, 6), sharex=True, sharey=True)
    for ax, m in zip(axes.flat, METHODS):
        if m in data:
            ax.plot(data[m][:, 1], data[m][:, 2], color=COLORS[m], lw=0.8)
        ax.plot(np.cos(theta), np.sin(theta), "k--", lw=0.7, alpha=0.45)
        ax.set_title(LABELS[m])
        ax.set_aspect("equal")
    fig.supxlabel("x")
    fig.supylabel("v")
    fig.suptitle("Phase portraits — harmonic oscillator  (dt = 0.05,  t = 50)", y=1.01)
    fig.tight_layout()
    save(fig, f"{out}/phase_portrait.png")

    # ── RK4 trajectory vs exact ───────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(9, 3.5))
    if "rk4" in data:
        t = data["rk4"][:, 0]
        ax.plot(t, data["rk4"][:, 1], color=COLORS["rk4"], label="RK4")
        ax.plot(t, np.cos(t), "k--", lw=0.9, alpha=0.7, label=r"$\cos(t)$ (exact)")
    ax.set_xlabel("t")
    ax.set_ylabel("x(t)")
    ax.set_title("RK4 trajectory vs exact solution")
    ax.legend()
    save(fig, f"{out}/trajectory.png")


# ══════════════════════════════════════════════════════════════════════════ #
# 2. Rössler attractor
# ══════════════════════════════════════════════════════════════════════════ #

def plot_rossler():
    out = "plots/rossler"
    os.makedirs(out, exist_ok=True)

    A, B, C = 0.2, 0.2, 5.7  # standard chaotic parameters

    # Long integration for the attractor (discard first 100 t for transient)
    full = run("rossler", A, B, C, 1, 0, 0, 0.01, 300, 1)
    if full.ndim != 2 or full.shape[1] < 4:
        return
    mask = full[:, 0] >= 100
    x, y, z = full[mask, 1], full[mask, 2], full[mask, 3]

    # ── Attractor projections ─────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axes[0].plot(x, y, lw=0.35, color="steelblue", alpha=0.85)
    axes[0].set_xlabel("x");  axes[0].set_ylabel("y")
    axes[0].set_title("x–y projection")
    axes[1].plot(x, z, lw=0.35, color="darkorange", alpha=0.85)
    axes[1].set_xlabel("x");  axes[1].set_ylabel("z")
    axes[1].set_title("x–z projection")
    fig.suptitle(f"Rössler attractor  (a={A}, b={B}, c={C})  —  RK4, dt = 0.01")
    fig.tight_layout()
    save(fig, f"{out}/attractor.png")

    # ── Time series ───────────────────────────────────────────────────────
    ts = run("rossler", A, B, C, 1, 0, 0, 0.01, 100, 1)
    fig, axes = plt.subplots(3, 1, figsize=(9, 7), sharex=True)
    for ax, col, lbl, clr in zip(axes, [1, 2, 3], ["x", "y", "z"],
                                  ["steelblue", "darkorange", "forestgreen"]):
        ax.plot(ts[:, 0], ts[:, col], color=clr, lw=0.9)
        ax.set_ylabel(lbl)
    axes[-1].set_xlabel("t")
    fig.suptitle(f"Rössler attractor — time series  (a={A}, b={B}, c={C},  RK4)")
    fig.tight_layout()
    save(fig, f"{out}/time_series.png")


# ══════════════════════════════════════════════════════════════════════════ #
# 3. Random walk
# ══════════════════════════════════════════════════════════════════════════ #

def plot_random_walk():
    out = "plots/random_walk"
    os.makedirs(out, exist_ok=True)

    # Unbiased walk: d=1, mu=0
    data = run("random_walk", 1000, 5000)
    if data.ndim != 2 or data.shape[1] < 4:
        return
    t, x_av, x2_av, x4_av = data[:, 0], data[:, 1], data[:, 2], data[:, 3]

    # ── Mean square displacement (log-log) ────────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 4.5))
    pos = t > 0
    ax.loglog(t[pos], x2_av[pos], color="steelblue", label=r"$\langle x^2 \rangle$")
    ax.loglog(t[pos], t[pos], "k--", lw=0.9, label=r"$\sim t$ (diffusive)")
    ax.set_xlabel("t")
    ax.set_ylabel(r"$\langle x^2 \rangle$")
    ax.set_title(r"Mean square displacement — 1D random walk  (d=1, $\mu$=0,  5000 traj.)")
    ax.legend()
    save(fig, f"{out}/mean_square_displacement.png")

    # ── Kurtosis → 3 (Gaussian limit) ────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 4.5))
    valid = t > 5
    kurt  = np.where(x2_av[valid] > 0,
                     x4_av[valid] / x2_av[valid] ** 2,
                     np.nan)
    ax.semilogx(t[valid], kurt, color="darkorange", label="kurtosis")
    ax.axhline(3.0, color="k", ls="--", lw=0.9, label="Gaussian limit (3)")
    ax.set_xlabel("t")
    ax.set_ylabel(r"$\langle x^4 \rangle \,/\, \langle x^2 \rangle^2$")
    ax.set_title("Kurtosis — convergence to Gaussian")
    ax.legend()
    save(fig, f"{out}/kurtosis.png")

    # ── Biased walk comparison ────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    mu_vals  = [0.0, 0.1, 0.3, 0.5]
    colors_m = ["steelblue", "darkorange", "forestgreen", "crimson"]
    for mu, c in zip(mu_vals, colors_m):
        d = run("random_walk", 500, 3000, 1.0 - mu**2, mu)
        if d.ndim != 2 or d.shape[1] < 3:
            continue
        t_m = d[:, 0]
        pos = t_m > 0
        axes[0].plot(t_m, d[:, 1], color=c, lw=1.1, label=rf"$\mu={mu}$")
        axes[1].loglog(t_m[pos], d[pos, 2], color=c, lw=1.1, label=rf"$\mu={mu}$")
    axes[0].set_xlabel("t");  axes[0].set_ylabel(r"$\langle x \rangle$")
    axes[0].set_title("Mean position (drift)")
    axes[0].legend()
    axes[1].set_xlabel("t");  axes[1].set_ylabel(r"$\langle x^2 \rangle$")
    axes[1].set_title("Mean square displacement")
    axes[1].legend()
    fig.suptitle("Effect of drift on 1D random walk  (3000 trajectories)")
    fig.tight_layout()
    save(fig, f"{out}/drift_comparison.png")


# ══════════════════════════════════════════════════════════════════════════ #
# 4. Percolation
# ══════════════════════════════════════════════════════════════════════════ #

def _perc_stats(L, n_cfg, rho_min, rho_max, d_rho):
    """Run percolation_ex and return (rhos, prob, mean_cluster_size) arrays."""
    data = run("percolation_ex", L, n_cfg, rho_min, rho_max, d_rho)
    if data.ndim != 2 or data.shape[1] < 5:
        return np.array([]), np.array([]), np.array([])
    rhos     = np.round(data[:, 1] / (L * L), 6)
    unique   = np.unique(rhos)
    probs, mcs = [], []
    for r in unique:
        mask = rhos == r
        probs.append(data[mask, 2].mean())
        mcs.append(data[mask, 4].mean())
    return unique, np.array(probs), np.array(mcs)


def plot_percolation():
    out = "plots/percolation"
    os.makedirs(out, exist_ok=True)
    print("    (running percolation — may take a few seconds)")

    specs  = [(10, 500, 0.02), (20, 200, 0.02), (40, 100, 0.02), (80, 50, 0.02)]
    colors = ["steelblue", "darkorange", "forestgreen", "crimson"]
    RHO_C  = 0.5927           # known critical density for square-lattice site percolation

    results = {}
    for (L, n_cfg, d_rho), c in zip(specs, colors):
        rhos, probs, mcs = _perc_stats(L, n_cfg, 0.3, 0.8, d_rho)
        results[L] = (rhos, probs, mcs, c)

    # ── Percolation probability ───────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for L, (rhos, probs, mcs, c) in results.items():
        ax.plot(rhos, probs, "o-", color=c, ms=3, lw=1.2, label=f"L = {L}")
    ax.axvline(RHO_C, color="gray", ls=":", lw=1.3,
               label=rf"$\rho_c \approx {RHO_C}$")
    ax.set_xlabel(r"$\rho$")
    ax.set_ylabel(r"$P(\,\mathrm{percolates}\,)$")
    ax.set_title("Site percolation probability — square lattice")
    ax.legend()
    save(fig, f"{out}/probability.png")

    # ── Mean cluster size (diverges at ρ_c) ──────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for L, (rhos, probs, mcs, c) in results.items():
        ax.semilogy(rhos, np.where(mcs > 0, mcs, np.nan),
                    "o-", color=c, ms=3, lw=1.2, label=f"L = {L}")
    ax.axvline(RHO_C, color="gray", ls=":", lw=1.3,
               label=rf"$\rho_c \approx {RHO_C}$")
    ax.set_xlabel(r"$\rho$")
    ax.set_ylabel(r"$\langle S \rangle$")
    ax.set_title("Mean cluster size — square lattice")
    ax.legend()
    save(fig, f"{out}/mean_cluster_size.png")


# ══════════════════════════════════════════════════════════════════════════ #
# 5. Lattice gas
# ══════════════════════════════════════════════════════════════════════════ #

def plot_lattice_gas():
    out = "plots/lattice_gas"
    os.makedirs(out, exist_ok=True)

    L     = 30
    steps = 2000
    rhos  = [0.05, 0.20, 0.40, 0.60]
    colors = ["steelblue", "darkorange", "forestgreen", "crimson"]

    fig, ax = plt.subplots(figsize=(8, 5))
    for rho, c in zip(rhos, colors):
        N    = int(rho * L * L)
        data = run("lattice_gas", L, N, steps)
        if data.ndim != 2 or data.shape[1] < 2:
            continue
        ax.plot(data[:, 0], data[:, 1], color=c, lw=1.2,
                label=rf"$\rho = {rho:.2f}$  (N = {N})")

    ax.axhline(1.0, color="gray", ls="--", lw=0.9, label="free diffusion (D = 1)")
    ax.set_xlabel("step")
    ax.set_ylabel(r"$\langle r^2 \rangle\,/\,\mathrm{step}$")
    ax.set_title(f"2D lattice gas — effective diffusion at different densities  (L = {L})")
    ax.legend()
    save(fig, f"{out}/diffusion.png")


# ══════════════════════════════════════════════════════════════════════════ #
# Entry point
# ══════════════════════════════════════════════════════════════════════════ #

if __name__ == "__main__":
    sections = [
        ("harmonic_oscillator", plot_harmonic_oscillator),
        ("rossler",             plot_rossler),
        ("random_walk",         plot_random_walk),
        ("percolation",         plot_percolation),
        ("lattice_gas",         plot_lattice_gas),
    ]
    for name, fn in sections:
        print(f"  {name}/")
        fn()
    print("Done.")
