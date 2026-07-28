"""Backward-compatible interface to the classical waiting-time simulation."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from QRSim.waiting_time import (
    accumulated_storage_rounds,
    completion_rounds,
    sample_continuum_waiting_times,
    sample_waiting_times,
)


class MonteCarloSim:
    """One realization of the repeater waiting-time process.

    New code should normally call :func:`sample_waiting_times`.  This class is
    retained for compatibility with the original QRSim scripts.
    """

    SUCCESS_PROBABILITY = 0.01
    DATA_DIRECTORY = Path(__file__).resolve().parent / "datas"
    CONTINUUM_SAMPLE_COUNT = 500_000
    CONTINUUM_SEED_BASE = 20_260_728
    CONTINUUM_DATA_VERSION = 1

    def __init__(
        self,
        link_num: int,
        cutoff_coefficient: float | int | None = None,
        p: float = SUCCESS_PROBABILITY,
        *,
        seed: int | None = None,
    ) -> None:
        if link_num < 1:
            raise ValueError("link_num must be positive")
        if not 0.0 < p <= 1.0:
            raise ValueError("p must lie in (0, 1]")
        if cutoff_coefficient is not None:
            raise NotImplementedError(
                "The old cutoff simulator does not implement the manuscript "
                "protocol reliably; cutoff restoration is tracked separately."
            )

        rng = np.random.default_rng(seed)
        self.link_num = link_num
        self.success_probability = p
        self.cut_off_coefficient = cutoff_coefficient
        self.success_rounds = rng.geometric(p, size=link_num).astype(np.int64)
        self.final_success_round = int(completion_rounds(self.success_rounds))
        self.total_storage_round = int(accumulated_storage_rounds(self.success_rounds))

    @staticmethod
    def simulate_bernoulli_trials(p: float, *, seed: int | None = None) -> int:
        """Return a geometric first-success round on ``{1, 2, ...}``."""

        if not 0.0 < p <= 1.0:
            raise ValueError("p must lie in (0, 1]")
        return int(np.random.default_rng(seed).geometric(p))

    def calculate_final_success_round(self) -> int:
        return int(completion_rounds(self.success_rounds))

    def calculate_total_storage_round(self) -> int:
        return int(accumulated_storage_rounds(self.success_rounds))

    @staticmethod
    def gen_freq_dist(
        link_num: int,
        success_prob: float,
        cutoff_coefficient: float | int | None = None,
        rnd_num: int = 10_000,
        *,
        seed: int | None = None,
        scaled: bool = True,
    ) -> list[np.ndarray]:
        """Generate completion/storage samples.

        Historical QRSim data store the scaled variables ``p Z_N`` and
        ``p W_N``.  ``scaled=True`` preserves that on-disk convention.
        """

        if cutoff_coefficient is not None:
            raise NotImplementedError(
                "Cutoff sampling must be restored from an explicit protocol definition."
            )
        z, w = sample_waiting_times(
            success_prob,
            link_count=link_num,
            trials=rnd_num,
            seed=seed,
        )
        if scaled:
            return [z * success_prob, w * success_prob]
        return [z, w]

    @staticmethod
    def save_freq_dist(
        datas: list[np.ndarray],
        repeater_nodes: int,
        cutoff_coefficient: int | float | None,
        *,
        success_prob: float,
        seed: int | None = None,
        output_directory: str | Path | None = None,
    ) -> Path:
        """Save scaled samples together with their provenance metadata."""

        if repeater_nodes < 0:
            raise ValueError("repeater_nodes cannot be negative")
        directory = Path(output_directory) if output_directory else MonteCarloSim.DATA_DIRECTORY
        directory.mkdir(parents=True, exist_ok=True)
        path = directory / f"distribution_n{repeater_nodes}_y{cutoff_coefficient}.npz"
        z, w = datas
        np.savez(
            path,
            Z=np.asarray(z),
            W=np.asarray(w),
            success_probability=float(success_prob),
            seed=-1 if seed is None else int(seed),
            repeater_nodes=int(repeater_nodes),
            link_count=int(repeater_nodes + 1),
            scaled=True,
        )
        return path

    @staticmethod
    def gen_continuum_dist(
        repeater_nodes: int,
        *,
        sample_count: int = CONTINUUM_SAMPLE_COUNT,
        seed: int | None = None,
    ) -> list[np.ndarray]:
        """Generate strict continuum samples from unit-rate exponentials."""

        z_bar, w_bar = sample_continuum_waiting_times(
            repeater_nodes,
            trials=sample_count,
            seed=seed,
        )
        return [z_bar, w_bar]

    @staticmethod
    def save_continuum_dist(
        datas: list[np.ndarray],
        repeater_nodes: int,
        *,
        seed: int,
        output_directory: str | Path | None = None,
    ) -> Path:
        """Cache reproducible continuum samples together with provenance."""

        if repeater_nodes < 0:
            raise ValueError("repeater_nodes cannot be negative")
        directory = Path(output_directory) if output_directory else MonteCarloSim.DATA_DIRECTORY
        directory.mkdir(parents=True, exist_ok=True)
        path = directory / f"distribution_n{repeater_nodes}_yNone.npz"
        z_bar, w_bar = datas
        np.savez(
            path,
            Z=np.asarray(z_bar, dtype=float),
            W=np.asarray(w_bar, dtype=float),
            sample_model="unit_rate_exponential_continuum",
            data_version=MonteCarloSim.CONTINUUM_DATA_VERSION,
            seed=int(seed),
            repeater_nodes=int(repeater_nodes),
            link_count=int(repeater_nodes + 1),
            sample_count=int(len(z_bar)),
        )
        return path

    @staticmethod
    def load_freq_dist(
        repeater_nodes: int,
        cutoff_coefficient: int | float | None = None,
        *,
        data_directory: str | Path | None = None,
        generate_if_missing: bool = True,
    ) -> np.lib.npyio.NpzFile:
        """Load or deterministically regenerate strict continuum samples.

        Large derived arrays are deliberately omitted from the reviewer
        package. A missing or legacy finite-probability cache is replaced
        using unit-rate exponential variables and a documented fixed seed.
        """

        directory = Path(data_directory) if data_directory else MonteCarloSim.DATA_DIRECTORY
        path = directory / f"distribution_n{repeater_nodes}_y{cutoff_coefficient}.npz"
        valid_cache = False
        if path.exists() and cutoff_coefficient is None:
            try:
                with np.load(path) as data:
                    valid_cache = (
                        str(data["sample_model"]) == "unit_rate_exponential_continuum"
                        and int(data["data_version"]) == MonteCarloSim.CONTINUUM_DATA_VERSION
                        and int(data["repeater_nodes"]) == repeater_nodes
                        and int(data["link_count"]) == repeater_nodes + 1
                        and int(data["seed"]) == MonteCarloSim.CONTINUUM_SEED_BASE + repeater_nodes
                        and int(data["sample_count"]) == MonteCarloSim.CONTINUUM_SAMPLE_COUNT
                        and len(data["Z"]) == MonteCarloSim.CONTINUUM_SAMPLE_COUNT
                        and len(data["W"]) == MonteCarloSim.CONTINUUM_SAMPLE_COUNT
                    )
            except (KeyError, ValueError):
                valid_cache = False

        if not valid_cache:
            if not generate_if_missing or cutoff_coefficient is not None:
                raise FileNotFoundError(path)
            seed = MonteCarloSim.CONTINUUM_SEED_BASE + repeater_nodes
            samples = MonteCarloSim.gen_continuum_dist(
                repeater_nodes,
                sample_count=MonteCarloSim.CONTINUUM_SAMPLE_COUNT,
                seed=seed,
            )
            MonteCarloSim.save_continuum_dist(
                samples,
                repeater_nodes,
                seed=seed,
                output_directory=directory,
            )
        return np.load(path)

    @staticmethod
    def plot_freq_dist(datas, link_num: int, cutoff_coefficient):
        """Plot cached continuum distributions (optional convenience helper)."""

        from matplotlib import pyplot as plt
        import seaborn as sns

        z = datas["Z"]
        w = datas["W"]
        sns.set(style="whitegrid")
        _, axes = plt.subplots(2, 1, figsize=(10, 6))
        sns.histplot(w, ax=axes[1], kde=False, binwidth=0.1, stat="density", color="skyblue")
        sns.histplot(z, ax=axes[0], kde=False, binwidth=0.1, stat="density", color="skyblue")
        axes[0].set(xlim=(0, 10), ylim=(0, 0.5), xlabel=r"$\overline{z}$")
        axes[1].set(xlim=(0, 10), ylim=(0, 0.5), xlabel=r"$\overline{w}$")
        axes[1].set_ylabel(r"$P\{\mathcal{W}_N\}$")
        plt.suptitle(f"Waiting-time distributions (links={link_num}, cutoff={cutoff_coefficient})")
        plt.show()


if __name__ == "__main__":
    data = MonteCarloSim.load_freq_dist(repeater_nodes=1)
    MonteCarloSim.plot_freq_dist(data, link_num=2, cutoff_coefficient=None)
