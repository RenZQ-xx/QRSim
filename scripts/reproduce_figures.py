"""Reproduce the manuscript's surviving computational figures.

Outputs are written to ``results/figures`` by default and never overwrite the
historical PNG files in the repository root.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

from QRSim.analytical import homogeneous_chain
from QRSim.MonteCarloSim import MonteCarloSim
from QRSim import HeraldedLink, RelayChain


MANUSCRIPT_PARAMETERS = {
    "detector_efficiency": 0.9,
    "dark_count_probability": 1.5e-6,
    "transmission_without_fiber": 0.08,
    "attenuation_db_per_km": 0.2,
    "visibility": 0.9,
    "coherence_time": 600.0,
    "swapping_quality": 0.995,
}

WAITING_TIME_SAMPLE_COUNT = MonteCarloSim.CONTINUUM_SAMPLE_COUNT
WAITING_TIME_SEED_BASE = MonteCarloSim.CONTINUUM_SEED_BASE


def ensure_waiting_time_data(repeater_nodes: set[int]) -> None:
    """Create the deterministic Monte Carlo samples needed by selected figures.

    Node count ``N`` uses seed ``20260728 + N`` and 500,000 draws from
    independent unit-rate exponential variables. Existing strict-continuum
    caches are retained; legacy finite-probability caches are regenerated.
    """

    for nodes in sorted(repeater_nodes):
        expected_seed = WAITING_TIME_SEED_BASE + nodes
        data = MonteCarloSim.load_freq_dist(repeater_nodes=nodes)
        if (
            str(data["sample_model"]) != "unit_rate_exponential_continuum"
            or int(data["seed"]) != expected_seed
            or len(data["W"]) != WAITING_TIME_SAMPLE_COUNT
        ):
            raise ValueError(f"invalid continuum cache for N={nodes}")


def save_figure(fig: plt.Figure, output: Path, name: str) -> Path:
    output.mkdir(parents=True, exist_ok=True)
    path = output / name
    fig.savefig(path, dpi=200, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(path.resolve())
    return path


def save_rows(output: Path, name: str, rows: list[dict]) -> Path:
    output.mkdir(parents=True, exist_ok=True)
    path = output / name
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def reproduce_fig2(output: Path) -> Path:
    """Distance scaling for N=1,3,5 with p_e=0.03."""

    distances = np.linspace(100, 1000, 100)
    colors = {1: "black", 3: "blue", 5: "green"}
    fig, axes = plt.subplots(1, 2, figsize=(11.36, 4.0))
    rows = []
    for nodes in (1, 3, 5):
        results = [
            homogeneous_chain(
                repeater_nodes=nodes,
                total_distance_km=float(distance),
                excitation_probability=0.03,
                **MANUSCRIPT_PARAMETERS,
            )
            for distance in distances
        ]
        axes[0].plot(distances, [r.fidelity for r in results], ".-", ms=3, lw=1, color=colors[nodes])
        axes[1].plot(
            distances,
            [r.rate for r in results],
            ".-",
            ms=3,
            lw=1,
            color=colors[nodes],
            label=f"N = {nodes}",
        )
        rows.extend(
            {
                "repeater_nodes": nodes,
                "total_distance_km": float(distance),
                "excitation_probability": 0.03,
                "fidelity": result.fidelity,
                "rate_hz": result.rate,
                "heralding_probability": result.elementary_link.heralding_probability,
                "coherence_factor": result.coherence_factor,
                "attempt_rate_hz": result.repetition_rate,
            }
            for distance, result in zip(distances, results)
        )
    axes[0].set(xlabel="Total distance (km)", ylabel="Fidelity", xlim=(100, 1000), ylim=(0, 1))
    axes[1].set(xlabel="Total distance (km)", ylabel="Rate (Hz)", xlim=(100, 1000), yscale="log")
    axes[1].yaxis.set_label_position("right")
    axes[1].yaxis.tick_right()
    axes[1].legend()
    fig.tight_layout()
    save_rows(output, "fig2_reproduced.csv", rows)
    return save_figure(fig, output, "fig2_reproduced.png")


def reproduce_fig3(output: Path) -> Path:
    """Fidelity--rate tradeoff at 1000 km for N=3,...,8."""

    # Use one common physical scan for every node count.  Curves below the
    # displayed fidelity window are clipped by the axes, as in the old figure.
    excitation_probabilities = np.linspace(0.001, 0.12, 120)
    colors = ["blue", "deepskyblue", "green", "gold", "red", "darkred"]
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    rows = []
    for nodes, color in zip(range(3, 9), colors):
        results = [
            homogeneous_chain(
                repeater_nodes=nodes,
                total_distance_km=1000,
                excitation_probability=float(pe),
                **MANUSCRIPT_PARAMETERS,
            )
            for pe in excitation_probabilities
        ]
        ax.plot(
            [r.rate for r in results],
            [r.fidelity for r in results],
            ".-",
            ms=3,
            lw=1,
            color=color,
            label=f"N = {nodes}",
        )
        rows.extend(
            {
                "repeater_nodes": nodes,
                "total_distance_km": 1000.0,
                "excitation_probability": float(pe),
                "fidelity": result.fidelity,
                "rate_hz": result.rate,
                "heralding_probability": result.elementary_link.heralding_probability,
                "coherence_factor": result.coherence_factor,
                "attempt_rate_hz": result.repetition_rate,
            }
            for pe, result in zip(excitation_probabilities, results)
        )
    ax.set(xlabel="Rate (Hz)", ylabel="Fidelity", xscale="log", xlim=(8e-3, 1), ylim=(0.45, 0.61))
    ax.legend()
    fig.tight_layout()
    save_rows(output, "fig3_reproduced.csv", rows)
    return save_figure(fig, output, "fig3_reproduced.png")


def reproduce_fig4_legacy(output: Path) -> Path:
    """Recover the numerical panels of historical fig4.

    This intentionally uses ``eta_0=0.10`` and ``T_coh=600 s`` because those
    are the values in the surviving script.  They conflict with the current
    caption/Sec. IV and are therefore labeled as a legacy reproduction.
    """

    cases = [
        (300.0, 0.024, 80, "blue"),
        (400.0, 0.041, 40, "green"),
        (500.0, 0.070, 30, "red"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    for total_distance, central_pe, point_count, color in cases:
        asymmetries = np.arange(point_count) * 0.5 / 80
        fidelities, rates, heralding = [], [], []
        for asymmetry in asymmetries:
            left_distance = total_distance * (0.5 - asymmetry)
            right_distance = total_distance - left_distance
            left_pe = central_pe * 10 ** (-0.2 * asymmetry * total_distance / 20)
            right_pe = min(1.0, central_pe * 10 ** (0.2 * asymmetry * total_distance / 20))
            links = []
            for distance, pe in ((left_distance, left_pe), (right_distance, right_pe)):
                link = HeraldedLink(
                    detector_efficiency=0.9,
                    dark_count_prob=1.5e-6,
                    transmission_without_fiber=0.10,
                    fiber_attenuation_coefficient=0.2,
                    distance=float(distance),
                    visibility=0.9,
                )
                link.generate_entanglement(float(pe))
                links.append(link)
            chain = RelayChain(Links=links, coherence_time=600)
            chain.solve(swapping_quality=0.995)
            fidelities.append(chain.final_fide)
            rates.append(chain.final_rate)
            heralding.append(np.mean([link.herald_prob for link in links]))
        label = rf"$L_{{tot}}={total_distance:.0f}$ km, $p_{{hr}}={heralding[0]:.1e}$"
        axes[0].plot(asymmetries, fidelities, color=color, label=label)
        axes[1].plot(asymmetries, rates, color=color, label=label)
    axes[0].set(xlabel=r"$|1/2-L_0/L_{tot}|$", ylabel="Fidelity", ylim=(0.2, 0.9))
    axes[1].set(xlabel=r"$|1/2-L_0/L_{tot}|$", ylabel="Rate (Hz)")
    axes[1].legend()
    fig.tight_layout()
    return save_figure(fig, output, "fig4_legacy_reproduced_panels.png")


def _solve_central_excitation_probability(total_distance: float, target_heralding: float) -> float:
    """Solve for p_e using the active Sec. IV parameter set."""

    low, high = 0.0, 1.0
    for _ in range(80):
        middle = (low + high) / 2
        result = homogeneous_chain(
            repeater_nodes=1,
            total_distance_km=total_distance,
            excitation_probability=middle,
            **{**MANUSCRIPT_PARAMETERS, "coherence_time": math.inf},
        )
        if result.elementary_link.heralding_probability < target_heralding:
            low = middle
        else:
            high = middle
    return (low + high) / 2


def _draw_fig4_schematic(axis: plt.Axes) -> None:
    """Draw the asymmetric one-repeater geometry used in fig4(a)."""

    node_style = dict(fill=False, linewidth=2.2)
    for x, color, label in ((0.05, "#0072B2", "Alice"), (0.48, "#7A3DB8", "Repeater"), (0.91, "#0072B2", "Bob")):
        axis.add_patch(Rectangle((x - 0.07, 0.58), 0.14, 0.22, edgecolor=color, **node_style))
        axis.text(x, 0.85, label, ha="center", va="bottom", fontsize=9)
    stations = ((0.27, r"$H_0$"), (0.73, r"$H_1$"))
    for x, label in stations:
        axis.add_patch(Rectangle((x - 0.045, 0.18), 0.09, 0.15, facecolor="#D9EAD3", edgecolor="#6AA84F", lw=1.8))
        axis.text(x, 0.255, label, ha="center", va="center", fontsize=9)
        axis.plot([x, x], [0.33, 0.69], "--", color="0.45", lw=1.5)
    axis.plot([0.12, 0.41], [0.69, 0.69], "--", color="0.45", lw=1.8)
    axis.plot([0.55, 0.84], [0.69, 0.69], "--", color="0.45", lw=1.8)
    axis.text(0.265, 0.73, r"$L_0$", ha="center", fontsize=10)
    axis.text(0.695, 0.73, r"$L_1$", ha="center", fontsize=10)
    axis.annotate("", xy=(0.48, 0.48), xytext=(0.40, 0.48), arrowprops=dict(arrowstyle="<->", lw=1.2))
    axis.text(0.44, 0.43, r"$\Delta L$", ha="center", fontsize=9)
    axis.set(xlim=(-0.05, 1.02), ylim=(0.08, 0.95))
    axis.set_aspect("equal")
    axis.axis("off")
    axis.set_title("(a)", loc="left")


def reproduce_fig4(output: Path) -> Path:
    """Asymmetric-placement panels using the active manuscript parameters."""

    cases = [
        (300.0, 1.4e-4, "blue"),
        (400.0, 7.7e-5, "green"),
        (500.0, 4.3e-5, "red"),
    ]
    fig = plt.figure(figsize=(15, 4))
    grid = fig.add_gridspec(1, 3, width_ratios=(1.15, 2, 2))
    schematic_axis = fig.add_subplot(grid[0, 0])
    axes = (fig.add_subplot(grid[0, 1]), fig.add_subplot(grid[0, 2]))
    _draw_fig4_schematic(schematic_axis)
    rows = []
    for total_distance, target_heralding, color in cases:
        central_pe = _solve_central_excitation_probability(total_distance, target_heralding)
        # The long-link memory has p_e = p_e(0) 10^(alpha b L_tot/20).
        # Stop at the physical boundary p_e=1 (and at b=1/2).
        max_asymmetry = min(
            0.5,
            math.log10(1 / central_pe) * 20 / (0.2 * total_distance),
        )
        asymmetries = np.linspace(0, max_asymmetry, 100)
        fidelities, rates = [], []
        for asymmetry in asymmetries:
            left_distance = total_distance * (0.5 - asymmetry)
            right_distance = total_distance - left_distance
            left_pe = central_pe * 10 ** (-0.2 * asymmetry * total_distance / 20)
            right_pe = min(1.0, central_pe * 10 ** (0.2 * asymmetry * total_distance / 20))
            links = []
            for distance, pe in ((left_distance, left_pe), (right_distance, right_pe)):
                link = HeraldedLink(
                    detector_efficiency=MANUSCRIPT_PARAMETERS["detector_efficiency"],
                    dark_count_prob=MANUSCRIPT_PARAMETERS["dark_count_probability"],
                    transmission_without_fiber=MANUSCRIPT_PARAMETERS["transmission_without_fiber"],
                    fiber_attenuation_coefficient=MANUSCRIPT_PARAMETERS["attenuation_db_per_km"],
                    distance=float(distance),
                    visibility=MANUSCRIPT_PARAMETERS["visibility"],
                )
                link.generate_entanglement(float(pe))
                links.append(link)
            chain = RelayChain(Links=links, coherence_time=math.inf)
            chain.solve(swapping_quality=MANUSCRIPT_PARAMETERS["swapping_quality"])
            fidelities.append(chain.final_fide)
            rates.append(chain.final_rate)
            rows.append(
                {
                    "total_distance_km": total_distance,
                    "target_central_heralding_probability": target_heralding,
                    "asymmetry": float(asymmetry),
                    "left_link_distance_km": float(left_distance),
                    "right_link_distance_km": float(right_distance),
                    "left_excitation_probability": float(left_pe),
                    "right_excitation_probability": float(right_pe),
                    "left_heralding_probability": float(links[0].herald_prob.real),
                    "right_heralding_probability": float(links[1].herald_prob.real),
                    "fidelity": float(chain.final_fide.real),
                    "rate_hz": float(chain.final_rate),
                }
            )
        label = rf"$L_{{tot}}={total_distance:.0f}$ km, $p_{{hr}}={target_heralding:.1e}$"
        axes[0].plot(asymmetries, fidelities, color=color, label=label)
        axes[1].plot(asymmetries, rates, color=color, label=label)
        print(
            f"fig4 L={total_distance:.0f} km: central p_e={central_pe:.8f}, "
            f"maximum asymmetry={max_asymmetry:.5f}"
        )
    axes[0].set(xlabel=r"$|1/2-L_0/L_{tot}|$", ylabel="Fidelity", xlim=(0, 0.5), ylim=(0.2, 0.9))
    axes[1].set(xlabel=r"$|1/2-L_0/L_{tot}|$", ylabel="Rate (Hz)", xlim=(0, 0.5))
    axes[0].set_title("(b)", loc="left")
    axes[1].set_title("(c)", loc="left")
    axes[1].legend()
    fig.tight_layout()
    save_rows(output, "fig4_reproduced.csv", rows)
    save_figure(fig, output, "fig4_reproduced.png")
    # Keep the historical filename as a compatibility pointer for earlier runs.
    return output / "fig4_reproduced.png"


def gamma_density(nodes: int, values: np.ndarray) -> np.ndarray:
    return (nodes + 1) * (1 - np.exp(-values)) ** nodes * np.exp(-values)


def pi_density(nodes: int, values: np.ndarray) -> np.ndarray:
    if nodes == 1:
        return np.exp(-values / 2) / 2
    if nodes == 2:
        return np.exp(-values / 4) / 9 + 2 * np.exp(-values / 2) / 3 - 7 * np.exp(-values) / 9
    if nodes == 3:
        return (
            4 * np.exp(-values / 4) / 15
            + (3 / 16 + values / 16) * np.exp(-values / 2)
            - 7 * np.exp(-values) / 6
            + 57 * np.exp(-3 * values / 2) / 80
        )
    raise ValueError("the manuscript gives explicit Pi densities only for N=1,2,3")


def reproduce_figure_xx(output: Path) -> Path:
    """Monte Carlo check of the analytical continuum densities."""

    grid = np.linspace(0, 10, 1000)
    bins = np.linspace(0, 10, 201)
    fig, axes = plt.subplots(2, 3, figsize=(10, 6), sharex=True, sharey=True)
    for column, nodes in enumerate((1, 2, 3)):
        data = MonteCarloSim.load_freq_dist(repeater_nodes=nodes)
        axes[0, column].hist(data["W"], bins=bins, density=True, color="skyblue", edgecolor="black", lw=0.25)
        axes[0, column].plot(grid, pi_density(nodes, grid), color="red", lw=1.5)
        axes[1, column].hist(data["Z"], bins=bins, density=True, color="skyblue", edgecolor="black", lw=0.25)
        axes[1, column].plot(grid, gamma_density(nodes, grid), color="red", lw=1.5)
        axes[0, column].set_title(f"N={nodes}")
        axes[0, column].set_xlim(0, 10)
        axes[0, column].set_ylim(0, 0.5)
        axes[0, column].set_xlabel(r"$\overline{w}$")
        axes[1, column].set_xlabel(r"$\overline{z}$")
    axes[0, 0].set_ylabel(r"$\Pi_N(\overline{w})$")
    axes[1, 0].set_ylabel(r"$\Gamma_N(\overline{z})$")
    fig.legend(["analytical"], loc="upper right")
    fig.tight_layout()
    return save_figure(fig, output, "Figure_xx_reproduced.png")


FIGURES = {
    "fig4": reproduce_fig4,
    "fig2": reproduce_fig2,
    "fig3": reproduce_fig3,
    "fig4_legacy": reproduce_fig4_legacy,
    "Figure_xx": reproduce_figure_xx,
}
DEFAULT_FIGURES = ("fig4", "fig2", "fig3", "Figure_xx")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("figures", nargs="*", metavar="FIGURE")
    parser.add_argument("--output", type=Path, default=Path("results/figures"))
    args = parser.parse_args()
    selected_figures = args.figures or DEFAULT_FIGURES
    unknown_figures = sorted(set(selected_figures) - FIGURES.keys())
    if unknown_figures:
        parser.error(
            "unknown figure name(s): "
            + ", ".join(unknown_figures)
            + "; choose from "
            + ", ".join(FIGURES)
        )
    required_waiting_time_data: set[int] = set()
    if "fig2" in selected_figures:
        required_waiting_time_data.update((1, 3, 5))
    if "fig3" in selected_figures:
        required_waiting_time_data.update(range(3, 9))
    if "Figure_xx" in selected_figures:
        required_waiting_time_data.update((1, 2, 3))
    ensure_waiting_time_data(required_waiting_time_data)
    for name in selected_figures:
        FIGURES[name](args.output)


if __name__ == "__main__":
    main()
