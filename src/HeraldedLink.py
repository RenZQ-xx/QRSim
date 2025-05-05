from cmath import sqrt, exp

from qutip import basis, tensor, ket2dm, Qobj, qeye

ket0, ket1 = [basis(2, i) for i in range(2)]
d00 = ket2dm(tensor(ket0, ket0))
d01 = ket2dm(tensor(ket0, ket1))
d10 = ket2dm(tensor(ket1, ket0))
d11 = ket2dm(tensor(ket1, ket1))
d10_01 = tensor(ket1, ket0) * tensor(ket0, ket1).dag()
d01_10 = tensor(ket0, ket1) * tensor(ket1, ket0).dag()
d00_11 = tensor(ket0, ket0) * tensor(ket1, ket1).dag()
d11_00 = tensor(ket1, ket1) * tensor(ket0, ket0).dag()
psi_up = (tensor(basis(2, 0), basis(2, 1)) + tensor(basis(2, 1), basis(2, 0))).unit()
psi_down = (tensor(basis(2, 0), basis(2, 1)) - tensor(basis(2, 1), basis(2, 0))).unit()
phi_up = (tensor(basis(2, 0), basis(2, 0)) + tensor(basis(2, 1), basis(2, 1))).unit()
phi_down = (tensor(basis(2, 0), basis(2, 0)) - tensor(basis(2, 1), basis(2, 1))).unit()
bell_pvm_set = {'psi+': psi_up, 'psi-': psi_down, 'phi+': phi_up, 'phi-': phi_down}


class HeraldedLink:
    LIGHT_SPEED: float = 200000.
    FIBER_ATTENUATION_COEFFICIENT: float = 0.2  # dB/km

    _next_link = None

    def __init__(self,
                 detector_efficiency: float or [float],
                 dark_count_prob: float or [float],
                 transmission_without_fiber: float or [float],  # excluding detector efficiency
                 distance: float or [float],
                 repetition_rate: float or None = None,
                 fiber_attenuation_coefficient: float or [float] = FIBER_ATTENUATION_COEFFICIENT,
                 visibility: float or [float] = 1.,
                 phase_shift: float = 0.,  # (rad)
                 scheme: str = "SPI",
                 excitation_prob: None | float or [float] = None,
                 coherence_time: None | float or [float] = None,
                 state_leakage_coefficient: float = 0.
                 ) -> None:

        self.state_leakage_coefficient = state_leakage_coefficient
        self.herald_state_assemble = None
        self.herald_prob_set = None
        self.initial_fidelity_set = None
        self.herald_state_set = None
        self.herald_prob = None
        self._scheme = scheme
        self._excitation_prob = [excitation_prob, excitation_prob] if type(
            excitation_prob) == float else excitation_prob
        self.detector_efficiency = [detector_efficiency, detector_efficiency] if type(
            detector_efficiency) == float else detector_efficiency
        self.dark_count_prob = [dark_count_prob, dark_count_prob] if type(dark_count_prob) == float else dark_count_prob
        transmission_without_fiber = [transmission_without_fiber, transmission_without_fiber] \
            if type(transmission_without_fiber) == float else transmission_without_fiber
        fiber_attenuation_coefficient = [fiber_attenuation_coefficient, fiber_attenuation_coefficient] \
            if type(fiber_attenuation_coefficient) == float else fiber_attenuation_coefficient
        distance = [distance / 2, distance / 2] if type(distance) == float else distance
        self.transmission = [
            transmission_without_fiber[0] * 10 ** (-fiber_attenuation_coefficient[0] * distance[0] / 10),
            transmission_without_fiber[0] * 10 ** (-fiber_attenuation_coefficient[0] * distance[0] / 10)]
        coherence_time = [coherence_time, coherence_time] if type(coherence_time) == float else coherence_time
        self.coherence_time = [float('inf'), float('inf')] if not coherence_time else coherence_time
        self.visibility = visibility
        self.phase_shift = phase_shift
        self.repetition_rate = 1 / (
                2 * max(distance) / HeraldedLink.LIGHT_SPEED + 15e-6) if not repetition_rate else repetition_rate

    @property
    def scheme(self):
        return self._scheme

    @scheme.setter
    def scheme(self, value):
        raise AttributeError("Cannot modify read-only attribute")

    @property
    def excitation_prob(self):
        return self._excitation_prob

    @excitation_prob.setter
    def excitation_prob(self, value):
        raise AttributeError("Cannot modify read-only attribute")

    @property
    def next_link(self):
        return self._next_link

    @next_link.setter
    def next_link(self, value):
        raise AttributeError("Cannot modify read-only attribute")

    def generate_entanglement(self, excitation_prob: None | float or [float] = None) -> bool:
        if type(excitation_prob) == list and len(excitation_prob) != 2:
            raise AttributeError("Invalid input parameters!")

        excitation_prob = [excitation_prob, excitation_prob] if type(excitation_prob) == float else excitation_prob
        if self.scheme == 'SPI':
            if excitation_prob[0] < 0 or excitation_prob[0] > 1 or excitation_prob[1] < 0 or excitation_prob[1] > 1:
                raise ValueError("invalid excitation probability")
            self._excitation_prob = excitation_prob
            self.generate_SPI_link()
        if self.scheme == 'TPI':
            if excitation_prob:
                raise ValueError("invalid input for TPI scheme")
            self.generate_TPI_link()
        self.herald_state_set["psi+"] = HeraldedLink.deal_with_phase_shift(self.herald_state_set["psi+"], self.phase_shift)
        self.herald_state_set["psi-"] = HeraldedLink.deal_with_phase_shift(self.herald_state_set["psi-"],
                                                                           self.phase_shift)

        self.deal_with_decoherence()
        self.initial_fidelity_set = {
            "psi+": HeraldedLink.calculate_fidelity(state=self.herald_state_set["psi+"], target=ket2dm(psi_up)),
            "psi-": HeraldedLink.calculate_fidelity(state=self.herald_state_set["psi-"], target=ket2dm(psi_down))}
        return True

    @staticmethod
    def deal_with_phase_shift(rho: Qobj, phase_shift: float) -> Qobj:
        return HeraldedLink.modify_non_diagonal_elements(rho, exp(complex(0, 1) * phase_shift))

    def deal_with_decoherence(self) -> None:
        t = 1 / self.repetition_rate
        self.herald_state_set["psi+"] = HeraldedLink.modify_non_diagonal_elements(self.herald_state_set["psi+"],
                                                                   exp(-t / self.coherence_time[0]) * exp(
                                                                       -t / self.coherence_time[1]))
        self.herald_state_set["psi-"] = HeraldedLink.modify_non_diagonal_elements(self.herald_state_set["psi-"],
                                                                     exp(-t / self.coherence_time[0]) * exp(
                                                                         -t / self.coherence_time[1]))

    @staticmethod
    def calculate_fidelity(state: Qobj, target: Qobj) -> float:
        if state.tr() == 0.:
            return 0
        F = (state * target * target.dag()).tr().real / state.tr().real
        return F

    # private function, please don't import
    def generate_SPI_link(self) -> bool:
        eta = self.transmission

        Pd = self.dark_count_prob
        eta_d = self.detector_efficiency
        V = self.visibility

        eta1, eta2 = eta[0], eta[1]
        pc1, pc2 = self.excitation_prob[0], self.excitation_prob[1]
        gamma1, gamma2 = pc1 * eta1, pc2 * eta2
        gamma = sqrt(gamma1 * gamma2).real
        Pr_true_0 = Pd
        Pr_true_1 = [1 - (1 - Pd[0]) * (1 - eta_d[0] * (1 - Pd[0])), 1 - (1 - Pd[1]) * (1 - eta_d[1] * (1 - Pd[1]))]
        Pr_true_2 = [1 - (1 - Pd[0]) * (1 - eta_d[0] * (1 - Pd[0])) ** 2,
                     1 - (1 - Pd[1]) * (1 - eta_d[1] * (1 - Pd[1])) ** 2]
        Pr_false_0 = [1 - Pd[0], 1 - Pd[1]]
        Pr_false_1 = [(1 - Pd[0]) * (1 - eta_d[0] * (1 - Pd[0])), (1 - Pd[1]) * (1 - eta_d[1] * (1 - Pd[1]))]
        Pr_false_2 = [(1 - Pd[0]) * (1 - eta_d[0] * (1 - Pd[0])) ** 2, (1 - Pd[1]) * (1 - eta_d[1] * (1 - Pd[1])) ** 2]

        F_01 = [Pr_false_0[0] * Pr_true_0[1],  # (0,0)
                Pr_false_0[0] * Pr_true_1[1],  # (0,1)
                Pr_false_0[0] * Pr_true_1[1],  # (0,1n)
                Pr_false_0[0] * Pr_true_2[1],  # (0,2)
                Pr_false_0[0] * Pr_true_2[1],  # (0,2n)
                Pr_false_1[0] * Pr_true_0[1],  # (1,0)
                Pr_false_1[0] * Pr_true_0[1],  # (1n,0)
                Pr_false_2[0] * Pr_true_0[1],  # (2,0)
                Pr_false_2[0] * Pr_true_0[1],  # (2n,0)
                Pr_false_1[0] * Pr_true_1[1],  # (1,1n)
                Pr_false_1[0] * Pr_true_1[1]]  # (1n,1)

        F_10 = [Pr_true_0[0] * Pr_false_0[1],  # (0,0)
                Pr_true_0[0] * Pr_false_1[1],  # (0,1)
                Pr_true_0[0] * Pr_false_1[1],  # (0,1n)
                Pr_true_0[0] * Pr_false_2[1],  # (0,2)
                Pr_true_0[0] * Pr_false_2[1],  # (0,2n)
                Pr_true_1[0] * Pr_false_0[1],  # (1,0)
                Pr_true_1[0] * Pr_false_0[1],  # (1n,0)
                Pr_true_2[0] * Pr_false_0[1],  # (2,0)
                Pr_true_2[0] * Pr_false_0[1],  # (2n,0)
                Pr_true_1[0] * Pr_false_1[1],  # (1,1n)
                Pr_true_1[0] * Pr_false_1[1]]  # (1n,1)

        diag: list[Qobj] = [(1 - pc1) * (1 - pc2) * d00 + (1 - eta2) * (1 - pc1) * pc2 * d01 \
                            + (1 - eta1) * pc1 * (1 - pc2) * d10 + (1 - eta1) * pc1 * (1 - eta2) * pc2 * d11,  # (0,0)
                            V / 2 * (1 - pc1) * gamma2 * d01 + V / 2 * gamma2 * (1 - eta1) * pc1 * d11 \
                            - sqrt(V) / 2 * sqrt((1 - pc1) * (1 - pc2)) * gamma * (d10_01 + d01_10) \
                            + 1 / 2 * gamma1 * (1 - pc2) * d10 + 1 / 2 * (gamma1 * (1 - eta2) * pc2) * d11,  # (0,1)
                            (1 - V) / 2 * ((1 - pc1) * gamma2 * d01 + gamma2 * (1 - eta1) * pc1 * d11),  # (0,1n)
                            V / 2 * gamma ** 2 * d11,  # (0,2)
                            (1 - V) / 4 * gamma ** 2 * d11,  # (0,2n)
                            V / 2 * (1 - pc1) * gamma2 * d01 + V / 2 * gamma2 * (1 - eta1) * pc1 * d11 \
                            + sqrt(V) / 2 * sqrt((1 - pc1) * (1 - pc2)) * gamma * (d10_01 + d01_10) \
                            + 1 / 2 * gamma1 * (1 - pc2) * d10 + 1 / 2 * gamma1 * (1 - eta2) * pc2 * d11,  # (1,0)
                            (1 - V) / 2 * ((1 - pc1) * gamma2 * d01 + (gamma2 * (1 - eta1) * pc1 * d11)),  # (1n,0)
                            V / 2 * gamma ** 2 * d11,  # (2,0)
                            (1 - V) / 4 * gamma ** 2 * d11,  # (2n,0)
                            (1 - V) / 4 * gamma ** 2 * d11,  # (1,1n)
                            (1 - V) / 4 * gamma ** 2 * d11]  # (1n,1)

        rho_D1_click = sum(F * d for F, d in zip(F_10, diag))
        # Psi+
        rho_D2_click = sum(F * d for F, d in zip(F_01, diag))
        # Psi-
        self.herald_prob_set = [rho_D1_click.tr(), rho_D2_click.tr()]
        self.herald_prob = rho_D1_click.tr() + rho_D2_click.tr()
        self.herald_state_set = {"psi+": rho_D1_click.unit(), "psi-": rho_D2_click.unit()}
        self.herald_state_assemble = {"psi+": rho_D1_click / self.herald_prob, "psi-": rho_D2_click / self.herald_prob}
        return True

    # private function, please don't import
    def generate_TPI_link(self,
                          r1=None,
                          r2=None) -> bool:
        if self.dark_count_prob[0] != self.dark_count_prob[1] or self.detector_efficiency[0] != \
                self.detector_efficiency[1]:
            raise AttributeError("TPI demands consistent detector parameters")
        if r1 is None:
            r1 = [1 / sqrt(2), 1 / sqrt(2)]
        if r2 is None:
            r2 = [1 / sqrt(2), 1 / sqrt(2)]
        alpha = r1[0]
        beta = r1[1]
        gamma = r2[0]
        delta = r2[1]

        V = self.visibility
        eta = self.transmission
        Pd = self.dark_count_prob
        eta_d = self.detector_efficiency

        eta1 = eta[0]
        eta2 = eta[1]

        Pr_true_0 = Pd
        Pr_true_1 = 1 - (1 - Pd) * (1 - eta_d * (1 - Pd))
        Pr_true_2 = 1 - (1 - Pd) * (1 - eta_d * (1 - Pd)) ** 2
        Pr_false_0 = 1 - Pd
        Pr_false_1 = (1 - Pd) * (1 - eta_d * (1 - Pd))
        Pr_false_2 = (1 - Pd) * (1 - eta_d * (1 - Pd)) ** 2

        alpha_x = alpha.conjugate()
        beta_x = beta.conjugate()
        gamma_x = gamma.conjugate()
        delta_x = delta.conjugate()

        diag_0 = (1 - eta1) * (1 - eta2) * (
                beta * delta * alpha_x * delta_x * d00 + beta * gamma * beta_x * gamma_x * d01
                + alpha * delta * alpha_x * delta_x * d10 + alpha * gamma * alpha_x * gamma_x * d11)
        diag_1 = (eta1 * (1 - eta2) + (1 - eta1) * eta2) * (beta * delta * alpha_x * delta_x * d00
                                                            + beta * gamma * beta_x * gamma_x * d01
                                                            + alpha * delta * alpha_x * delta_x * d10
                                                            + alpha * gamma * alpha_x * gamma_x * d11)
        diag_2 = (1 + V) / 2 * eta1 * eta2 * (beta * delta * alpha_x * delta_x * d00
                                              + alpha * gamma * alpha_x * gamma_x * d11)
        diag_3 = 1 / 2 * eta1 * eta2 * (
                beta * gamma * beta_x * gamma_x * d01 + V * beta * gamma * alpha_x * delta_x * d01_10
                + V * alpha * delta * beta_x * gamma_x * d10_01 + alpha * delta * alpha_x * delta_x * d10)
        diag_4 = 1 / 2 * eta1 * eta2 * (
                beta * gamma * beta_x * gamma_x * d01 - V * beta * gamma * alpha_x * delta_x * d01_10
                - V * alpha * delta * beta_x * gamma_x * d10_01 + alpha * delta * alpha_x * delta_x * d10)
        diag_5 = (1 - V) / 2 * eta1 * eta2 * (beta * delta * alpha_x * delta_x * d00
                                              + alpha * gamma * alpha_x * gamma_x * d11)

        Prcl_0 = 2 * Pr_true_0 ** 2 * Pr_false_0 ** 2
        Prcl_1 = Pr_true_0 * Pr_true_1 * Pr_false_0 * Pr_false_0 + Pr_true_0 * Pr_true_0 * Pr_false_0 * Pr_false_1
        Prcl_2 = Pr_true_0 * Pr_true_2 * Pr_false_0 * Pr_false_0 + Pr_true_0 * Pr_true_0 * Pr_false_0 * Pr_false_2
        Prcl_3 = Pr_true_1 * Pr_true_1 * Pr_false_0 * Pr_false_0 + Pr_true_0 * Pr_true_0 * Pr_false_1 * Pr_false_1
        Prcl_4 = 2 * Pr_true_0 * Pr_true_1 * Pr_false_0 * Pr_false_1

        rho_psi_up = Prcl_0 * diag_0 \
                     + Prcl_1 * diag_1 \
                     + Prcl_2 * diag_2 \
                     + Prcl_3 * diag_3 \
                     + Prcl_4 * (diag_4 + diag_5)
        rho_psi_down = Prcl_0 * diag_0 \
                       + Prcl_1 * diag_1 \
                       + Prcl_2 * diag_2 \
                       + Prcl_3 * diag_4 \
                       + Prcl_4 * (diag_3 + diag_5)

        self.herald_prob = rho_psi_up.tr() + rho_psi_down.tr()
        self.herald_prob_set = [rho_psi_up.tr(), rho_psi_down.tr()]
        self.herald_state_set = {"psi+": rho_psi_up.unit(), "psi-": rho_psi_down.unit()}
        self.herald_state_assemble = {"psi+": rho_psi_up / self.herald_prob, "psi-": rho_psi_down / self.herald_prob}
        return True

    @staticmethod
    def modify_non_diagonal_elements(rho: Qobj, multiplier: float | complex) -> Qobj:
        if rho.dims != [[2, 2], [2, 2]]:
            raise AttributeError("Invalid input matrix! (dimension [[2, 2], [2, 2]] is demand)")
        p00 = (rho * d00).tr().real
        p01 = (rho * d01).tr().real
        p10 = (rho * d10).tr().real
        p11 = (rho * d11).tr().real
        p10_01 = (rho * d10_01).tr().real
        p01_10 = (rho * d01_10).tr().real
        p00_11 = (rho * d00_11).tr().real
        p11_00 = (rho * d11_00).tr().real

        p10_01 = p10_01 * multiplier.conjugate()
        p01_10 = p01_10 * multiplier
        p00_11 = p00_11 * multiplier
        p11_00 = p11_00 * multiplier.conjugate()
        rho = sum(P * d for P, d in zip([p00, p01, p10, p11, p10_01, p01_10, p00_11, p11_00],
                                        [d00, d01, d10, d11, d10_01, d01_10, d00_11, d11_00]))

        return rho




if __name__ == '__main__':
    
    link_12 = HeraldedLink(
        detector_efficiency=0.9,
        dark_count_prob=1.5e-6,
        transmission_without_fiber=0.08,
        distance=200.,
        visibility=0.9,
        phase_shift=0.05,
        scheme="SPI",
        coherence_time=10.)

    link_12.generate_entanglement(excitation_prob=0.03)
    print(link_12.herald_prob)
    print(link_12.herald_state_set)



    link_34 = HeraldedLink(
        detector_efficiency=[0.9, 0.8],
        dark_count_prob=[1.5e-6, 1.0e-6],
        transmission_without_fiber=0.08,
        distance=[110., 90.],
        visibility=0.9,
        phase_shift=0.05,
        scheme="SPI",
        coherence_time=[10., 8.])
    link_34.generate_entanglement(excitation_prob=[0.03, 0.02])
    print(link_34.herald_prob_set)
    print(link_34.herald_state_set)
    print(link_34.initial_fidelity_set)

