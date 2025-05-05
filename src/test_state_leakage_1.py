from RelayChain import *

from MonteCarloSim import *

gamma = 0.06
alpha = 0.07
alpha1 = alpha * (1 - gamma) / (1 - alpha * gamma)
alpha2 = 0.
N = 6
L = 1000.


link_1 = HeraldedLink.HeraldedLink(
    detector_efficiency=1.,
    dark_count_prob=1.5e-6,
    transmission_without_fiber=0.08,
    distance=L / (N+1),
    visibility=0.93,
    phase_shift=0.05,
    scheme="SPI",
    coherence_time=600.)
link_1.generate_entanglement(excitation_prob=alpha1)
print(link_1.herald_state_assemble)
print(link_1.herald_prob)
link_2 = HeraldedLink.HeraldedLink(
    detector_efficiency=1.,
    dark_count_prob=1.5e-6,
    transmission_without_fiber=0.08,
    distance=L / (N+1),
    visibility=0.93,
    phase_shift=0.05,
    scheme="SPI",
    coherence_time=600.)
link_2.generate_entanglement(excitation_prob=alpha2)
print(link_2.initial_fidelity_set)
link_12 = HeraldedLink.HeraldedLink(
    detector_efficiency=1.,
    dark_count_prob=1.5e-6,
    transmission_without_fiber=0.08,
    distance=L / (N+1),
    visibility=0.93,
    phase_shift=0.05,
    scheme="SPI",
    coherence_time=600.)
link_12.generate_entanglement(excitation_prob=[alpha1,alpha2])
print(link_12.initial_fidelity_set)
print("----------------\n")
link_new = HeraldedLink.HeraldedLink(
    detector_efficiency=1.,
    dark_count_prob=1.5e-6,
    transmission_without_fiber=0.08,
    distance=L/ (N+1),
    visibility=0.93,
    phase_shift=0.05,
    scheme="SPI",
    coherence_time=600.)
link_new.herald_prob = (1-alpha*gamma) **2 * link_1.herald_prob + (alpha*gamma)**2 * link_2.herald_prob \
                      + 2*(1-alpha*gamma)*alpha*gamma * link_12.herald_prob
print(link_new.herald_prob)
link_new.herald_state_assemble = {keys: (1-alpha*gamma) **2 * link_1.herald_state_assemble[keys] + (alpha*gamma)**2 * link_2.herald_state_assemble[keys] \
                      + 2*(1-alpha*gamma)*alpha*gamma * link_12.herald_state_assemble[keys] for keys in ['psi+','psi-']}
print(link_new.herald_state_assemble)

chain12 = RelayChain(Link=link_new, N=N)
chain12.solve(swapping_quality=0.99)
print(chain12.final_fidelity_without_decoherence)
print(chain12.final_rate)

