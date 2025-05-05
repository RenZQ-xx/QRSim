from qutip import Qobj, qeye, tensor, qzero

import HeraldedLink
from HeraldedLink import bell_pvm_set
from src.MonteCarloSim import MonteCarloSim


class RelayChain:
    zero_state = tensor(qzero(2), qeye(2)) * 0

    def __init__(self,
                 Link: HeraldedLink.HeraldedLink | None = None,
                 N: int = 0,
                 Links=None
                 ) -> None:
        if not RelayChain.isvalid(Link=Link, Links=Links):
            raise AttributeError("Invalid input attribute!")
        self.node_num = N if Link else len(Links) - 1
        self.scheme = Link.scheme if Link else Links[0].scheme
        import copy
        self.link_list = [copy.copy(Link) for _ in range(self.node_num + 1)] if Link else [copy.copy(Links[i]) for i in
                                                                                           range(self.node_num + 1)]

        self.initial_states = Link.herald_state_set if Link else [Links[i].herald_state_set for i in
                                                                  range(self.node_num + 1)]
        self.initial_fides = Link.initial_fidelity_set if Link else [Links[i].initial_fidelity_set for i in
                                                                     range(self.node_num + 1)]
        self.herald_prob = Link.herald_prob if Link else [Links[i].herald_prob for i in range(self.node_num + 1)]

        self.repetition_rate = min([link.repetition_rate for link in self.link_list])
        self.isSolved = False

        self.final_state = None
        self.final_fide = None
        self.final_rate = None

        self.final_state_without_decoherence = None
        self.final_fidelity_without_decoherence = None
        self.final_fidelity_set_without_decoherence = None

    def set_excitation_prob(self, excitation_prob: float or [float]) -> None:
        excitation_prob = [excitation_prob] * (self.node_num + 1) if type(excitation_prob) == float else excitation_prob
        if len(excitation_prob) != self.node_num + 1:
            raise AttributeError
        for i in range(self.node_num + 1):
            self.link_list[i].generate_entanglement(excitation_prob=excitation_prob[i])

    def solve(self,
              initial_Bell: None or str = None,
              BSM_results: None or str or [str] = None,
              swapping_quality: float or [float] = 1.):
        swapping_quality = [swapping_quality] * self.node_num if type(swapping_quality) == float else swapping_quality

        if not initial_Bell:
            initial_Bell = [self.link_list[i].herald_state_assemble for i in range(self.node_num + 1)]
        elif type(initial_Bell) == str:
            initial_Bell = [{'psi+': RelayChain.zero_state, 'psi-': RelayChain.zero_state,
                             initial_Bell: self.link_list[i].herald_state_assemble[initial_Bell]}
                            for i in range(self.node_num + 1)]

        if not BSM_results:
            BSM_results = [bell_pvm_set for _ in range(self.node_num)]
        elif type(BSM_results) == str:
            BSM_results = [{'psi+': RelayChain.zero_state,
                            'psi-': RelayChain.zero_state,
                            'phi+': RelayChain.zero_state,
                            'phi-': RelayChain.zero_state,
                            BSM_results: bell_pvm_set[BSM_results]} for _ in range(self.node_num)]
        else:
            BSM_results = [dict({'psi+': RelayChain.zero_state,
                                 'psi-': RelayChain.zero_state,
                                 'phi+': RelayChain.zero_state,
                                 'phi-': RelayChain.zero_state},
                                **{BSM_results[i]: bell_pvm_set[BSM_results[i]]}) for i in range(self.node_num)]

        self.final_state_without_decoherence = RelayChain.swap_n_links(initial_Bell, BSM_results, swapping_quality)
        self.final_fidelity_set_without_decoherence = {
            keys: HeraldedLink.HeraldedLink.calculate_fidelity(self.final_state_without_decoherence[keys],
                                                               bell_pvm_set[keys])
            for keys in self.final_state_without_decoherence.keys()}
        self.final_fidelity_without_decoherence = sum([self.final_state_without_decoherence[keys].tr()
                                                       * HeraldedLink.HeraldedLink.calculate_fidelity(self.final_state_without_decoherence[keys],
                                                                                                      bell_pvm_set[keys]) for keys in
                                                       self.final_state_without_decoherence.keys()])

        datas = MonteCarloSim.load_freq_dist(link_num=self.node_num + 1, cutoff_coefficient=None)
        Z = datas['Z']
        W = datas['W']
        EZ = sum(Z) / len(Z)

        self.final_rate = self.repetition_rate * self.herald_prob / EZ

        self.isSolved = True

    @staticmethod
    def endure_depolarising_noise(rho_in: Qobj, lmd: float) -> Qobj:
        rho_out = lmd * rho_in \
                  + (1 - lmd) / 16 * qeye([2, 2, 2, 2]) * rho_in.tr()
        return rho_out

    @staticmethod
    def swap_n_links(initial_Bell: [dict],
                     BSM_results: [dict],
                     swapping_quality: [float]) -> [dict]:
        n = len(initial_Bell)

        temp_state = dict({'phi+': RelayChain.zero_state, 'phi-': RelayChain.zero_state}, **initial_Bell[0])
        for i in range(1, n):
            temp_state = RelayChain.swap_2links(temp_state, initial_Bell[i], BSM_results[i - 1],
                                                swapping_quality[i - 1])
        return temp_state

    @staticmethod
    def swap_2links(initial_Bell_1: dict,
                    initial_Bell_2: dict,
                    BSM_results: dict,
                    swapping_quality: float):

        element1 = tensor(initial_Bell_1['psi+'], initial_Bell_2['psi+']) + tensor(initial_Bell_1['psi-'],
                                                                                   initial_Bell_2['psi-'])
        element2 = tensor(initial_Bell_1['psi+'], initial_Bell_2['psi-']) + tensor(initial_Bell_1['psi-'],
                                                                                   initial_Bell_2['psi+'])
        element3 = tensor(initial_Bell_1['phi+'], initial_Bell_2['psi+']) + tensor(initial_Bell_1['phi-'],
                                                                                   initial_Bell_2['psi-'])
        element4 = tensor(initial_Bell_1['phi+'], initial_Bell_2['psi-']) + tensor(initial_Bell_1['phi-'],
                                                                                   initial_Bell_2['psi+'])

        element1 = RelayChain.endure_depolarising_noise(element1, swapping_quality)
        element2 = RelayChain.endure_depolarising_noise(element2, swapping_quality)
        element3 = RelayChain.endure_depolarising_noise(element3, swapping_quality)
        element4 = RelayChain.endure_depolarising_noise(element4, swapping_quality)

        BSM_results = {keys: tensor(qeye(2), BSM_results[keys] * BSM_results[keys].dag(), qeye(2)) for keys in
                       BSM_results.keys()}
        swapped_Bell_psi_up = (element1 * BSM_results['psi+']).ptrace([0, 3]) \
                              + (element2 * BSM_results['psi-']).ptrace([0, 3]) \
                              + (element3 * BSM_results['phi+']).ptrace([0, 3]) \
                              + (element4 * BSM_results['phi-']).ptrace([0, 3])
        swapped_Bell_psi_down = (element2 * BSM_results['psi+']).ptrace([0, 3]) \
                                + (element1 * BSM_results['psi-']).ptrace([0, 3]) \
                                + (element4 * BSM_results['phi+']).ptrace([0, 3]) \
                                + (element3 * BSM_results['phi-']).ptrace([0, 3])
        swapped_Bell_phi_up = (element3 * BSM_results['psi+']).ptrace([0, 3]) \
                              + (element4 * BSM_results['psi-']).ptrace([0, 3]) \
                              + (element1 * BSM_results['phi+']).ptrace([0, 3]) \
                              + (element2 * BSM_results['phi-']).ptrace([0, 3])
        swapped_Bell_phi_down = (element4 * BSM_results['psi+']).ptrace([0, 3]) \
                                + (element3 * BSM_results['psi-']).ptrace([0, 3]) \
                                + (element2 * BSM_results['phi+']).ptrace([0, 3]) \
                                + (element1 * BSM_results['phi-']).ptrace([0, 3])

        return {'psi+': swapped_Bell_psi_up, 'psi-': swapped_Bell_psi_down,
                'phi+': swapped_Bell_phi_up, 'phi-': swapped_Bell_phi_down}

    @staticmethod
    def isvalid(Link: HeraldedLink.HeraldedLink or None = None,
                Links: [HeraldedLink.HeraldedLink] or None = None,
                ) -> bool:
        if (Link and Links) or not (Link or Links):
            return False
        elif Links:
            for i in range(1, len(Links)):
                if Links[i].scheme != Links[0].scheme:
                    return False
        return True


if __name__ == '__main__':
    from HeraldedLink import HeraldedLink

    link_12 = HeraldedLink(
        detector_efficiency=0.9,
        dark_count_prob=1.5e-6,
        transmission_without_fiber=0.08,
        distance=200.)

    link_12.generate_entanglement(0.03)

    chain = RelayChain(Link=link_12,
                       N=1,
                       Links=None)
    chain.solve(swapping_quality=0.99)
