import random
from matplotlib import pyplot as plt
import seaborn as sns

import numpy as np


class MonteCarloSim:
    SUCCESS_PROBABILITY = 0.01
    success_rounds = []
    n = 0
    y = 0

    def __init__(self,
                 link_num: int,
                 cutoff_coefficient: None | float | int,
                 p=SUCCESS_PROBABILITY):
        self.cut_off_coefficient = cutoff_coefficient
        self.link_num = link_num

        MonteCarloSim.init_class_val(link_num, p, cutoff_coefficient)
        self.final_success_round = MonteCarloSim.calculate_final_success_round()
        self.total_storage_round = MonteCarloSim.calculate_final_success_round()

    @classmethod
    def init_class_val(cls, link_num: int, p: float, cutoff_coefficient: float):
        cls.success_rounds = [0] * link_num
        cls.SUCCESS_PROBABILITY = p
        cls.n = link_num
        cls.y = cutoff_coefficient

    @staticmethod
    def simulate_bernoulli_trials(p):
        trials = 0
        while random.random() > p:
            trials += 1
        return trials

    @classmethod
    def do_monte_carlo_sim(cls, p=SUCCESS_PROBABILITY):
        n = cls.n
        y = cls.y
        if y != 0:
            cls.success_rounds = [MonteCarloSim.simulate_bernoulli_trials(p) for _ in range(n)]

        else:
            link_ages: [int] = [-1] * n
            times_of_cut_off = 0
            final_success_round = 0
            while not all([link_ages[i] >= 0 for i in range(n)]):
                for i in range(n):
                    if link_ages[i] < 0:
                        if random.random() < p:
                            link_ages[i] = 0
                    elif p * link_ages[i] > y:
                        link_ages[i] = -1
                        times_of_cut_off += 1
                    else:
                        link_ages[i] += 1
                final_success_round += 1
            cls.success_rounds = [final_success_round - link_ages[i] for i in range(n)]

    @classmethod
    def calculate_final_success_round(cls):
        return max(cls.success_rounds)

    @classmethod
    def calculate_total_storage_round(cls):
        n = cls.n
        index = cls.success_rounds.index(max(cls.success_rounds))

        if index == 0 or index == n - 1:
            return 2 * (max(cls.success_rounds) - min(cls.success_rounds))
        else:
            A1 = min(cls.success_rounds[:index])
            A2 = min(cls.success_rounds[index + 1:])
            return 2 * ((max(cls.success_rounds) - A1) + (max(cls.success_rounds) - A2))

    @staticmethod
    def gen_freq_dist(link_num: int,
                      success_prob: float,
                      cutoff_coefficient: None | float | int = None,
                      rnd_num: int = 10000):
        Z = []
        W = []
        for i in range(rnd_num):
            trial = MonteCarloSim(link_num=link_num,
                                  cutoff_coefficient=cutoff_coefficient,
                                  p=success_prob)
            Z.append(trial.final_success_round * success_prob)
            W.append(trial.final_success_round * success_prob)
        return [Z, W]

    @staticmethod
    def save_freq_dist(datas: [[float]], link_num: int, cutoff_coefficient: int | float | None):
        Z, W = datas[0], datas[1]
        np.savez('distribution_n{}_y{}.npz'.format(link_num, cutoff_coefficient), Z=Z, W=W)

    @staticmethod
    def load_freq_dist(link_num: int, cutoff_coefficient: int | float | None):
        datas = np.load('..\datas\distribution_n{}_y{}.npz'.format(link_num, cutoff_coefficient))
        return datas


    @staticmethod
    def plot_freq_dist(datas, link_num, cutoff_coefficient):
        Z = datas['Z']
        W = datas['W']
        sns.set(style="whitegrid")

        # 绘制直方图，设置适当的 bins
        fig, axes = plt.subplots(2, 1, figsize=(10, 6))
        sns.histplot(W, ax=axes[1], kde=False, binwidth=0.1, stat="density", color="skyblue", linewidth=0.5)
        sns.histplot(Z, ax=axes[0], kde=False, binwidth=0.1, stat="density", color="skyblue", linewidth=0.5)

        axes[0].set_xlim([0, 10])
        axes[0].set_xticks([0, 2, 4, 6, 8, 10])
        axes[0].set_ylim([0, 0.5])
        axes[0].set_xlabel('$\overline{z}$')

        axes[1].set_xlim([0, 10])
        axes[1].set_ylim([0, 0.5])
        axes[1].set_xticks([0, 2, 4, 6, 8, 10])
        axes[1].set_xlabel('$\overline{w}$')
        axes[1].set_ylabel('$P\{\mathcal{W}_N\}$')

        plt.title(f"Probability Distribution of Results (n={link_num}, y={cutoff_coefficient})", fontsize=16)
        # 显示图形
        plt.show()


if __name__ == "__main__":
    data = MonteCarloSim.load_freq_dist(link_num=1, cutoff_coefficient=None)
    MonteCarloSim.plot_freq_dist(data, link_num=1, cutoff_coefficient=None)

