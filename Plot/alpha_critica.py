import argparse
import csv
import glob
import os
import re
import shutil
import subprocess

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import ConnectionPatch, Rectangle


def step_number(path):
    match = re.search(r"_(\d+)\.dat$", os.path.basename(path))
    return int(match.group(1)) if match else -1


def read_last_bisection_pair(path):
    with open(path, "r", encoding="utf-8") as file:
        rows = list(csv.DictReader(file))
    if not rows:
        raise ValueError(f"No hay datos en {path}")

    last = rows[-1]
    return float(last["weak"]), float(last["strong"])


def run_case(exe, nr, rmax, tfinal, phi0, outdir):
    os.makedirs(outdir, exist_ok=True)
    exe_name = os.path.basename(exe)
    local_exe = os.path.join(outdir, exe_name)
    shutil.copy2(exe, local_exe)

    for path in glob.glob(os.path.join(outdir, "CEesferico1D_*.dat")):
        os.remove(path)

    command = [local_exe, str(nr), str(rmax), str(tfinal), f"{phi0:.17g}"]
    log_path = os.path.join(outdir, "run.log")
    with open(log_path, "w", encoding="utf-8") as log:
        subprocess.run(command, cwd=outdir, stdout=log, stderr=subprocess.STDOUT, check=True)


def load_snapshot_times(log_path):
    times = {}
    if not os.path.exists(log_path):
        return times

    pattern = re.compile(r"(?:Paso|Colapso geometrico en n =)\s+(\d+),\s+t =\s+([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)")
    with open(log_path, "r", encoding="utf-8") as file:
        for line in file:
            match = pattern.search(line)
            if match:
                times[int(match.group(1))] = parse_fortran_float(match.group(2))
    return times


def load_alpha_central(pattern, dt, step_parity, time_by_step=None):
    files = sorted(glob.glob(pattern), key=step_number)
    if not files:
        raise FileNotFoundError(f"No encontre archivos con el patron: {pattern}")

    time_by_step = time_by_step or {}
    times = []
    alpha = []
    for path in files:
        step = step_number(path)
        if step_parity is not None and step % 2 != step_parity:
            continue
        data = load_fortran_table(path)
        times.append(time_by_step.get(step, step*dt))
        alpha.append(data[0, 5])

    return np.array(times), np.array(alpha)


def stop_after_alpha_below(series, threshold):
    if threshold is None:
        return series

    times, alpha = series
    crossing = np.flatnonzero(alpha < threshold)
    if crossing.size == 0:
        return series

    stop = int(crossing[0]) + 1
    return times[:stop], alpha[:stop]


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


def auto_zoom_limits(weak, strong):
    t_w, a_w = weak
    t_s, a_s = strong
    rise_start = float(t_w[np.argmin(a_w)])

    def window(before, after):
        xmin = max(0.0, rise_start - before)
        xmax = rise_start + after
        mask_w = (t_w >= xmin) & (t_w <= xmax)
        mask_s = (t_s >= xmin) & (t_s <= xmax)
        values = np.concatenate([a_w[mask_w], a_s[mask_s]])
        ymin = float(np.min(values))
        ymax = float(np.max(values))
        margin = max(0.01, 0.10*(ymax - ymin))
        return [xmin, xmax, ymin - margin, 0.5]

    return window(0.45, 1.35), window(0.15, 0.80)


def plot_zoomed_pair(weak, strong, weak_phi, strong_phi, output, zoom1, zoom2):
    t_w, a_w = weak
    t_s, a_s = strong

    plt.rcParams.update(
        {
            "font.family": "serif",
            "mathtext.fontset": "cm",
            "axes.unicode_minus": False,
        }
    )

    if zoom1 is None or zoom2 is None:
        auto_zoom1, auto_zoom2 = auto_zoom_limits(weak, strong)
        zoom1 = auto_zoom1 if zoom1 is None else zoom1
        zoom2 = auto_zoom2 if zoom2 is None else zoom2

    fig = plt.figure(figsize=(9.2, 6.0))
    grid = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.95], hspace=0.32, wspace=0.28)
    ax_top = fig.add_subplot(grid[0, :])
    ax_left = fig.add_subplot(grid[1, 0])
    ax_right = fig.add_subplot(grid[1, 1])

    label_w = rf"$p$ = {weak_phi:.18f}"
    label_s = rf"$p$ = {strong_phi:.18f}"
    style_w = dict(color="blue", linewidth=1.4, label=label_w)
    style_s = dict(color="darkorange", linewidth=1.4, linestyle=(0, (4, 4)), label=label_s)

    for axis in (ax_top, ax_left, ax_right):
        axis.plot(t_w, a_w, **style_w)
        axis.plot(t_s, a_s, **style_s)
        axis.grid(True, color="black", linestyle=":", linewidth=0.7)
        axis.tick_params(direction="in", top=True, right=True, labelsize=15)

    ax_top.set_ylabel(r"$\alpha_\mathrm{central}$", fontsize=22)
    ax_top.legend(loc="lower right", fontsize=10, frameon=True, fancybox=False, framealpha=1.0)

    for axis, limits in ((ax_left, zoom1), (ax_right, zoom2)):
        xmin, xmax, ymin, ymax = limits
        axis.set_xlim(xmin, xmax)
        axis.set_ylim(ymin, ymax)

    ax_left.set_xlabel("")
    ax_right.set_xlabel("")
    fig.supxlabel(r"$t$", fontsize=22, y=0.03)

    for limits, target_axis in ((zoom1, ax_left), (zoom2, ax_right)):
        xmin, xmax, ymin, ymax = limits
        rect = Rectangle(
            (xmin, ymin),
            xmax - xmin,
            ymax - ymin,
            fill=False,
            edgecolor="0.72",
            linestyle=(0, (4, 6)),
            linewidth=1.0,
        )
        ax_top.add_patch(rect)
        for source_x, target_x in ((xmin, 0.0), (xmax, 1.0)):
            con = ConnectionPatch(
                xyA=(source_x, ymin),
                coordsA=ax_top.transData,
                xyB=(target_x, 1.0),
                coordsB=target_axis.transAxes,
                color="0.72",
                linestyle=(0, (4, 6)),
                linewidth=1.0,
            )
            fig.add_artist(con)

    fig.savefig(output, dpi=220, bbox_inches="tight")
    print(f"Grafica guardada en: {output}")


def parse_limits(text):
    values = [float(value) for value in text.split(",")]
    if len(values) != 4:
        raise argparse.ArgumentTypeError("Usa xmin,xmax,ymin,ymax")
    return values


def parse_args():
    parser = argparse.ArgumentParser(
        description="Grafica alpha_central(t) para el par weak/strong de biseccion.csv."
    )
    parser.add_argument("--csv", default="biseccion_phi0_fortran/biseccion.csv")
    parser.add_argument("--exe", default="CEesferico1D.exe")
    parser.add_argument("--nr", type=int, default=640)
    parser.add_argument("--rmax", type=float, default=64.0)
    parser.add_argument("--tfinal", type=float, default=20.0)
    parser.add_argument("--weak-dir", default="alpha_critica_weak")
    parser.add_argument("--strong-dir", default="alpha_critica_strong")
    parser.add_argument("--skip-run", action="store_true", help="No corre el exe; solo lee los .dat existentes.")
    parser.add_argument("--output", default="alpha_central_critica.png")
    parser.add_argument("--zoom1", type=parse_limits, default=None)
    parser.add_argument("--zoom2", type=parse_limits, default=None)
    parser.add_argument(
        "--step-parity",
        type=int,
        choices=(0, 1),
        default=None,
        help="Grafica solo pasos pares (0) o impares (1). Util para quitar oscilaciones par/impar.",
    )
    parser.add_argument(
        "--stop-below",
        type=float,
        default=None,
        help="Corta cada curva justo despues de que alpha_central baja de este valor.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    weak_phi, strong_phi = read_last_bisection_pair(args.csv)
    dt = 0.5*args.rmax/args.nr

    if not args.skip_run:
        run_case(args.exe, args.nr, args.rmax, args.tfinal, weak_phi, args.weak_dir)
        run_case(args.exe, args.nr, args.rmax, args.tfinal, strong_phi, args.strong_dir)

    weak_times = load_snapshot_times(os.path.join(args.weak_dir, "run.log"))
    strong_times = load_snapshot_times(os.path.join(args.strong_dir, "run.log"))
    weak = load_alpha_central(os.path.join(args.weak_dir, "CEesferico1D_*.dat"), dt, args.step_parity, weak_times)
    strong = load_alpha_central(os.path.join(args.strong_dir, "CEesferico1D_*.dat"), dt, args.step_parity, strong_times)
    weak = stop_after_alpha_below(weak, args.stop_below)
    strong = stop_after_alpha_below(strong, args.stop_below)
    plot_zoomed_pair(weak, strong, weak_phi, strong_phi, args.output, args.zoom1, args.zoom2)


if __name__ == "__main__":
    main()
