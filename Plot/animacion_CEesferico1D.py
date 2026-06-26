import argparse
import glob
import os
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation


COLUMNS = {
    "scalar": 1,
    "Phi": 2,
    "Pi": 3,
    "a": 4,
    "alpha": 5,
}


def step_number(path):
    match = re.search(r"_(\d+)\.dat$", os.path.basename(path))
    return int(match.group(1)) if match else -1


def load_snapshots(pattern, variable):
    files = sorted(glob.glob(pattern), key=step_number)
    if not files:
        raise FileNotFoundError(f"No encontre archivos con el patron: {pattern}")

    column = COLUMNS[variable]
    snapshots = []

    for path in files:
        data = load_fortran_table(path)
        snapshots.append(
            {
                "file": path,
                "step": step_number(path),
                "r": data[:, 0],
                "y": data[:, column],
            }
        )

    return snapshots


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


def make_animation(snapshots, variable, interval, output, xmax, ymin, ymax):
    reference = max(snapshots, key=lambda snapshot: len(snapshot["r"]))
    r = reference["r"]
    y_frames = []

    for snapshot in snapshots:
        if len(snapshot["r"]) == len(r) and np.allclose(snapshot["r"], r):
            y_frames.append(snapshot["y"])
        else:
            y_frames.append(np.interp(r, snapshot["r"], snapshot["y"]))

    y_all = np.array(y_frames)

    if xmax is not None:
        mask = r <= xmax
        r = r[mask]
        y_all = y_all[:, mask]

    y_min = float(np.min(y_all)) if ymin is None else ymin
    y_max = float(np.max(y_all)) if ymax is None else ymax
    margin = 0.08*(y_max - y_min) if y_max > y_min else 1.0

    if variable == "a":
        if xmax is None:
            xmax = 16.0
            mask = r <= xmax
            r = r[mask]
            y_all = y_all[:, mask]
        if ymin is None:
            ymin = 0.9
        if ymax is None:
            ymax = 2.0
        y_min = ymin
        y_max = ymax
        margin = 0.0

    fig_size = (8.6, 2.8) if variable == "a" else (8, 4.8)
    line_color = "tab:orange" if variable == "a" else None

    fig, ax = plt.subplots(figsize=fig_size)
    line, = ax.plot(r, y_all[0], color=line_color, lw=1.8)
    title = ax.set_title("")

    xlabel = r"$r$" if variable == "a" else "r"
    ylabel = r"$a(t,r)$" if variable == "a" else variable

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(float(np.min(r)), float(np.max(r)))
    ax.set_ylim(y_min - margin, y_max + margin)
    if variable == "a":
        title.set_visible(False)
        ax.grid(True, color="black", linestyle=":", linewidth=0.6, alpha=0.7)
        ax.tick_params(direction="in", top=True, right=True)
        fig.tight_layout(pad=0.4)
    else:
        ax.grid(True, alpha=0.3)

    def update(frame):
        snapshot = snapshots[frame]
        line.set_ydata(y_all[frame])
        if variable != "a":
            title.set_text(f"{variable}(r), paso n = {snapshot['step']}")
        return line, title

    animation = FuncAnimation(
        fig,
        update,
        frames=len(snapshots),
        interval=interval,
        blit=False,
        repeat=True,
    )

    if output:
        extension = os.path.splitext(output)[1].lower()
        if extension == ".gif":
            animation.save(output, writer="pillow", fps=max(1, 1000//interval))
        elif extension == ".mp4":
            animation.save(output, writer="ffmpeg", fps=max(1, 1000//interval))
        else:
            raise ValueError("La salida debe terminar en .gif o .mp4")
        print(f"Animacion guardada en: {output}")
    else:
        plt.show()


def save_or_show_figure(fig, output):
    if output:
        extension = os.path.splitext(output)[1].lower()
        if extension not in [".png", ".pdf", ".svg"]:
            raise ValueError("Para --mode time la salida debe terminar en .png, .pdf o .svg")
        fig.savefig(output, dpi=180, bbox_inches="tight")
        print(f"Grafica guardada en: {output}")
    else:
        fig.savefig("grafica_tiempo.png", dpi=180, bbox_inches="tight")
        print("Grafica guardada en: grafica_tiempo.png")


def make_panel_plot(snapshots, variable, dt, output, xmax, ymin, ymax, ncols):
    if dt is None:
        dt = 1.0

    nframes = len(snapshots)
    ncols = max(1, ncols)
    nrows = int(np.ceil(nframes/ncols))

    label = r"\alpha" if variable == "alpha" else variable

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(3.0*ncols, 2.15*nrows),
        sharex=True,
        sharey=True,
        squeeze=False,
    )

    for axis in axes.flat:
        axis.set_visible(False)

    for axis, snapshot in zip(axes.flat, snapshots):
        r = snapshot["r"]
        y = snapshot["y"]

        if xmax is not None:
            mask = r <= xmax
            r = r[mask]
            y = y[mask]

        axis.set_visible(True)
        axis.plot(r, y, color="blue", lw=1.0)
        axis.set_title(rf"$t = {snapshot['step']*dt:.2f}$", y=0.84, fontsize=13)
        axis.grid(True, color="black", linestyle=":", linewidth=0.6)
        axis.tick_params(direction="in", top=True, right=True, labelsize=10)

        if xmax is not None:
            axis.set_xlim(0.0, xmax)
        if ymin is not None or ymax is not None:
            axis.set_ylim(ymin, ymax)

    for row in range(nrows):
        axes[row, 0].set_ylabel(rf"${label}(t,r)$", fontsize=12)
    for col in range(ncols):
        axes[nrows - 1, col].set_xlabel(r"$r$", fontsize=12)

    if output:
        extension = os.path.splitext(output)[1].lower()
        if extension not in [".png", ".pdf", ".svg"]:
            raise ValueError("Para --mode panels la salida debe terminar en .png, .pdf o .svg")
        fig.savefig(output, dpi=220, bbox_inches="tight")
        print(f"Paneles guardados en: {output}")
    else:
        fig.savefig("paneles_alpha.png", dpi=220, bbox_inches="tight")
        print("Paneles guardados en: paneles_alpha.png")


def make_time_plot(snapshots, variable, r0, dt, output, ymin, ymax, tmax):
    times, values = time_series_at_radius(snapshots, r0, dt)
    times, values = truncate_time_series(times, values, tmax)

    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.plot(times, values, marker="o", lw=1.8, ms=3)

    xlabel = "t" if dt is not None else "n"
    ax.set_xlabel(xlabel)
    ax.set_ylabel(f"{variable}(r={r0:g})")
    ax.set_title(f"Evolucion temporal de {variable} en r = {r0:g}")
    ax.grid(True, alpha=0.3)

    if ymin is not None or ymax is not None:
        ax.set_ylim(ymin, ymax)

    save_or_show_figure(fig, output)


def time_series_at_radius(snapshots, r0, dt):
    times = []
    values = []

    for snapshot in snapshots:
        r = snapshot["r"]
        y = snapshot["y"]
        value = float(np.interp(r0, r, y))
        time = snapshot["step"]*dt if dt is not None else snapshot["step"]
        times.append(time)
        values.append(value)

    return np.array(times, dtype=float), np.array(values, dtype=float)


def truncate_time_series(times, values, tmax):
    if tmax is None:
        return times, values

    mask = times <= tmax
    if not np.any(mask):
        raise ValueError(f"No hay datos con tiempo menor o igual a {tmax}.")
    return times[mask], values[mask]


def make_time_animation(snapshots, variable, r0, dt, interval, output, ymin, ymax, tmax):
    times, values = time_series_at_radius(snapshots, r0, dt)
    times, values = truncate_time_series(times, values, tmax)

    fig, ax = plt.subplots(figsize=(8, 4.8))
    line, = ax.plot([], [], lw=2.0)
    point, = ax.plot([], [], marker="o", ms=5)
    title = ax.set_title("")

    xlabel = "t" if dt is not None else "n"
    ax.set_xlabel(xlabel)
    ax.set_ylabel(f"{variable}(r={r0:g})")
    ax.set_xlim(float(np.min(times)), float(np.max(times)))

    y_min = float(np.min(values)) if ymin is None else ymin
    y_max = float(np.max(values)) if ymax is None else ymax
    margin = 0.08*(y_max - y_min) if y_max > y_min else 1.0
    ax.set_ylim(y_min - margin, y_max + margin)
    ax.grid(True, alpha=0.3)

    def update(frame):
        line.set_data(times[:frame + 1], values[:frame + 1])
        point.set_data([times[frame]], [values[frame]])
        title.set_text(
            f"Evolucion temporal de {variable} en r = {r0:g}, "
            f"{xlabel} = {times[frame]:.4g}"
        )
        return line, point, title

    animation = FuncAnimation(
        fig,
        update,
        frames=len(times),
        interval=interval,
        blit=False,
        repeat=True,
    )

    if output:
        extension = os.path.splitext(output)[1].lower()
        if extension == ".gif":
            animation.save(output, writer="pillow", fps=max(1, 1000//interval))
        elif extension == ".mp4":
            animation.save(output, writer="ffmpeg", fps=max(1, 1000//interval))
        else:
            raise ValueError("Para --mode time-animation la salida debe terminar en .gif o .mp4")
        print(f"Animacion guardada en: {output}")
    else:
        animation.save("animacion_tiempo.gif", writer="pillow", fps=max(1, 1000//interval))
        print("Animacion guardada en: animacion_tiempo.gif")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Anima las salidas radiales de CEesferico1D."
    )
    parser.add_argument(
        "--pattern",
        default="CEesferico1D_*.dat",
        help="Patron de archivos .dat a leer.",
    )
    parser.add_argument(
        "--variable",
        choices=sorted(COLUMNS),
        default="scalar",
        help="Variable a animar.",
    )
    parser.add_argument(
        "--interval",
        type=int,
        default=250,
        help="Tiempo entre cuadros en milisegundos.",
    )
    parser.add_argument(
        "--mode",
        choices=["profile", "time", "time-animation", "panels"],
        default="profile",
        help="profile anima variable(r); time grafica variable(t); time-animation anima variable(t); panels hace subgraficas.",
    )
    parser.add_argument(
        "--r0",
        type=float,
        default=0.0,
        help="Radio donde evaluar la variable cuando --mode time.",
    )
    parser.add_argument(
        "--dt",
        type=float,
        default=None,
        help="Paso temporal fisico. Si se omite en --mode time, usa n en el eje horizontal.",
    )
    parser.add_argument(
        "--tmax",
        type=float,
        default=None,
        help="Tiempo maximo a mostrar en --mode time o --mode time-animation.",
    )
    parser.add_argument(
        "--output",
        default="",
        help="Archivo de salida opcional, por ejemplo animacion.gif o animacion.mp4.",
    )
    parser.add_argument(
        "--xmax",
        type=float,
        default=None,
        help="Valor maximo de r para hacer zoom en la region central.",
    )
    parser.add_argument(
        "--ymin",
        type=float,
        default=None,
        help="Limite inferior opcional del eje vertical.",
    )
    parser.add_argument(
        "--ymax",
        type=float,
        default=None,
        help="Limite superior opcional del eje vertical.",
    )
    parser.add_argument(
        "--every",
        type=int,
        default=1,
        help="Usar una de cada N salidas para --mode panels.",
    )
    parser.add_argument(
        "--max-panels",
        type=int,
        default=9,
        help="Numero maximo de paneles en --mode panels.",
    )
    parser.add_argument(
        "--ncols",
        type=int,
        default=3,
        help="Numero de columnas en --mode panels.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    snapshots = load_snapshots(args.pattern, args.variable)
    if args.mode == "time":
        make_time_plot(
            snapshots,
            args.variable,
            args.r0,
            args.dt,
            args.output,
            args.ymin,
            args.ymax,
            args.tmax,
        )
    elif args.mode == "time-animation":
        make_time_animation(
            snapshots,
            args.variable,
            args.r0,
            args.dt,
            args.interval,
            args.output,
            args.ymin,
            args.ymax,
            args.tmax,
        )
    elif args.mode == "panels":
        selected = snapshots[::max(1, args.every)][:args.max_panels]
        make_panel_plot(
            selected,
            args.variable,
            args.dt,
            args.output,
            args.xmax,
            args.ymin,
            args.ymax,
            args.ncols,
        )
    else:
        make_animation(
            snapshots,
            args.variable,
            args.interval,
            args.output,
            args.xmax,
            args.ymin,
            args.ymax,
        )


if __name__ == "__main__":
    main()
