import argparse
import glob
import os
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation


def step_number(path):
    match = re.search(r"_(\d+)\.dat$", os.path.basename(path))
    return int(match.group(1)) if match else -1


def parse_fortran_float(token):
    try:
        return float(token)
    except ValueError:
        fixed = re.sub(r"(\d)([+-]\d+)$", r"\1E\2", token)
        return float(fixed)


def load_fortran_table(path):
    rows = []
    with open(path, "r", encoding="utf-8") as file:
        for line in file:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            rows.append([parse_fortran_float(token) for token in line.split()])
    return np.array(rows, dtype=float)


def load_mass_snapshots(pattern, every):
    files = sorted(glob.glob(pattern), key=step_number)
    if not files:
        raise FileNotFoundError(f"No encontre archivos con el patron: {pattern}")

    snapshots = []
    for path in files[::max(1, every)]:
        data = load_fortran_table(path)
        snapshots.append(
            {
                "file": path,
                "step": step_number(path),
                "r": data[:, 0],
                "mass": data[:, 1],
                "compactness": data[:, 2],
            }
        )

    return snapshots


def interpolate_to_reference_grid(snapshots, key):
    reference = max(snapshots, key=lambda snapshot: len(snapshot["r"]))
    r = reference["r"]
    frames = []

    for snapshot in snapshots:
        if len(snapshot["r"]) == len(r) and np.allclose(snapshot["r"], r):
            frames.append(snapshot[key])
        else:
            frames.append(np.interp(r, snapshot["r"], snapshot[key]))

    return r, np.array(frames)


def apply_xmax(r, *arrays, xmax):
    if xmax is None:
        return (r, *arrays)

    mask = r <= xmax
    return (r[mask], *(array[:, mask] for array in arrays))


def y_limits(values, ymin, ymax):
    y_min = float(np.min(values)) if ymin is None else ymin
    y_max = float(np.max(values)) if ymax is None else ymax
    margin = 0.08*(y_max - y_min) if y_max > y_min else 1.0
    return y_min - margin, y_max + margin


def make_animation(snapshots, interval, output, xmax, mass_ymin, mass_ymax, comp_ymin, comp_ymax):
    r, mass_frames = interpolate_to_reference_grid(snapshots, "mass")
    _, comp_frames = interpolate_to_reference_grid(snapshots, "compactness")
    r, mass_frames, comp_frames = apply_xmax(r, mass_frames, comp_frames, xmax=xmax)

    plt.rcParams.update(
        {
            "font.family": "serif",
            "mathtext.fontset": "cm",
            "axes.unicode_minus": False,
        }
    )

    fig, (ax_mass, ax_comp) = plt.subplots(
        2,
        1,
        figsize=(8.0, 6.2),
        sharex=True,
        constrained_layout=True,
    )

    mass_line, = ax_mass.plot(r, mass_frames[0], color="blue", linewidth=1.6)
    comp_line, = ax_comp.plot(r, comp_frames[0], color="darkorange", linewidth=1.6)
    title = fig.suptitle("", fontsize=14)

    ax_mass.set_ylabel(r"$m_\mathrm{MS}(r)$", fontsize=16)
    ax_comp.set_ylabel(r"$2m_\mathrm{MS}/r$", fontsize=16)
    ax_comp.set_xlabel(r"$r$", fontsize=16)

    ax_mass.set_ylim(*y_limits(mass_frames, mass_ymin, mass_ymax))
    ax_comp.set_ylim(*y_limits(comp_frames, comp_ymin, comp_ymax))
    ax_comp.set_xlim(float(np.min(r)), float(np.max(r)))

    for axis in (ax_mass, ax_comp):
        axis.grid(True, color="black", linestyle=":", linewidth=0.7)
        axis.tick_params(direction="in", top=True, right=True, labelsize=12)

    def update(frame):
        snapshot = snapshots[frame]
        mass_line.set_ydata(mass_frames[frame])
        comp_line.set_ydata(comp_frames[frame])
        title.set_text(f"Masa de Misner-Sharp, paso n = {snapshot['step']}")
        return mass_line, comp_line, title

    animation = FuncAnimation(
        fig,
        update,
        frames=len(snapshots),
        interval=interval,
        blit=False,
        repeat=True,
    )

    extension = os.path.splitext(output)[1].lower()
    if extension == ".gif":
        animation.save(output, writer="pillow", fps=max(1, 1000//interval))
    elif extension == ".mp4":
        animation.save(output, writer="ffmpeg", fps=max(1, 1000//interval))
    else:
        raise ValueError("La salida debe terminar en .gif o .mp4")

    print(f"Animacion guardada en: {output}")


def make_panel_plot(snapshots, output, xmax):
    r, mass_frames = interpolate_to_reference_grid(snapshots, "mass")
    _, comp_frames = interpolate_to_reference_grid(snapshots, "compactness")
    r, mass_frames, comp_frames = apply_xmax(r, mass_frames, comp_frames, xmax=xmax)

    nframes = len(snapshots)
    fig, axes = plt.subplots(
        nframes,
        2,
        figsize=(8.0, 2.0*nframes),
        sharex=True,
        squeeze=False,
        constrained_layout=True,
    )

    for i, snapshot in enumerate(snapshots):
        axes[i, 0].plot(r, mass_frames[i], color="blue", linewidth=1.1)
        axes[i, 1].plot(r, comp_frames[i], color="darkorange", linewidth=1.1)
        axes[i, 0].set_ylabel(rf"$n={snapshot['step']}$", fontsize=11)

    axes[0, 0].set_title(r"$m_\mathrm{MS}(r)$")
    axes[0, 1].set_title(r"$2m_\mathrm{MS}/r$")
    axes[-1, 0].set_xlabel(r"$r$")
    axes[-1, 1].set_xlabel(r"$r$")

    for axis in axes.flat:
        axis.grid(True, color="black", linestyle=":", linewidth=0.6)
        axis.tick_params(direction="in", top=True, right=True, labelsize=10)

    extension = os.path.splitext(output)[1].lower()
    if extension not in [".png", ".pdf", ".svg"]:
        raise ValueError("Para --mode panels la salida debe terminar en .png, .pdf o .svg")
    fig.savefig(output, dpi=220, bbox_inches="tight")
    print(f"Paneles guardados en: {output}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Anima o grafica m_MS(r) y 2m/r desde archivos mass_*.dat."
    )
    parser.add_argument("--pattern", default="mass_*.dat")
    parser.add_argument("--output", default="animacion_masa.gif")
    parser.add_argument("--interval", type=int, default=250)
    parser.add_argument("--xmax", type=float, default=None)
    parser.add_argument("--every", type=int, default=1, help="Usar una de cada N salidas.")
    parser.add_argument(
        "--mode",
        choices=["animation", "panels"],
        default="animation",
        help="animation crea gif/mp4; panels crea una figura estatica.",
    )
    parser.add_argument("--max-panels", type=int, default=8)
    parser.add_argument("--mass-ymin", type=float, default=None)
    parser.add_argument("--mass-ymax", type=float, default=None)
    parser.add_argument("--comp-ymin", type=float, default=None)
    parser.add_argument("--comp-ymax", type=float, default=None)
    return parser.parse_args()


def main():
    args = parse_args()
    snapshots = load_mass_snapshots(args.pattern, args.every)

    if args.mode == "panels":
        make_panel_plot(snapshots[:args.max_panels], args.output, args.xmax)
    else:
        make_animation(
            snapshots,
            args.interval,
            args.output,
            args.xmax,
            args.mass_ymin,
            args.mass_ymax,
            args.comp_ymin,
            args.comp_ymax,
        )


if __name__ == "__main__":
    main()
