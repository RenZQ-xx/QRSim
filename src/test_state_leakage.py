from RelayChain import *
from HeraldedLink import HeraldedLink
from MonteCarloSim import *

# \gamma = 0.
# 200km
link_12 = HeraldedLink(
        detector_efficiency=1.,
        dark_count_prob=1.5e-6,
        transmission_without_fiber=0.08,
        distance=200./2,
        visibility=0.93,
        phase_shift=0.05,
        scheme="SPI",
        coherence_time=600.)
link_12.generate_entanglement(excitation_prob=0.03)
print(link_12.initial_fidelity_set)
chain = RelayChain(Link=link_12,
                   N=1)
chain.solve(swapping_quality=0.99)
print(chain.final_fidelity_without_decoherence)
print(chain.final_rate)

link_12 = HeraldedLink(
        detector_efficiency=1.,
        dark_count_prob=1.5e-6,
        transmission_without_fiber=0.08,
        distance=200./3,
        visibility=0.93,
        phase_shift=0.05,
        scheme="SPI",
        coherence_time=600.)
link_12.generate_entanglement(excitation_prob=0.03)
print(link_12.initial_fidelity_set)
chain = RelayChain(Link=link_12,
                   N=2)
chain.solve(swapping_quality=0.99)
print(chain.final_fidelity_without_decoherence)
print(chain.final_rate)

link_12 = HeraldedLink(
        detector_efficiency=1.,
        dark_count_prob=1.5e-6,
        transmission_without_fiber=0.08,
        distance=200./4,
        visibility=0.93,
        phase_shift=0.05,
        scheme="SPI",
        coherence_time=600.)
link_12.generate_entanglement(excitation_prob=0.03)
print(link_12.initial_fidelity_set)
chain = RelayChain(Link=link_12,
                   N=3)
chain.solve(swapping_quality=0.99)
print(chain.final_fidelity_without_decoherence)
print(chain.final_rate)

# 500km
link_12 = HeraldedLink(
        detector_efficiency=1.,
        dark_count_prob=1.5e-6,
        transmission_without_fiber=0.08,
        distance=500./2,
        visibility=0.93,
        phase_shift=0.05,
        scheme="SPI",
        coherence_time=600.)
link_12.generate_entanglement(excitation_prob=0.03)
print(link_12.initial_fidelity_set)
chain = RelayChain(Link=link_12,
                   N=1)
chain.solve(swapping_quality=0.99)
print(chain.final_fidelity_without_decoherence)
print(chain.final_rate)

link_12 = HeraldedLink(
        detector_efficiency=1.,
        dark_count_prob=1.5e-6,
        transmission_without_fiber=0.08,
        distance=500./3,
        visibility=0.93,
        phase_shift=0.05,
        scheme="SPI",
        coherence_time=600.)
link_12.generate_entanglement(excitation_prob=0.03)
print(link_12.initial_fidelity_set)
chain = RelayChain(Link=link_12,
                   N=2)
chain.solve(swapping_quality=0.99)
print(chain.final_fidelity_without_decoherence)
print(chain.final_rate)

link_12 = HeraldedLink(
        detector_efficiency=1.,
        dark_count_prob=1.5e-6,
        transmission_without_fiber=0.08,
        distance=500./4,
        visibility=0.93,
        phase_shift=0.05,
        scheme="SPI",
        coherence_time=600.)
link_12.generate_entanglement(excitation_prob=0.03)
print(link_12.initial_fidelity_set)
chain = RelayChain(Link=link_12,
                   N=3)
chain.solve(swapping_quality=0.99)
print(chain.final_fidelity_without_decoherence)
print(chain.final_rate)

# 800km
link_12 = HeraldedLink(
        detector_efficiency=1.,
        dark_count_prob=1.5e-6,
        transmission_without_fiber=0.08,
        distance=800./2,
        visibility=0.93,
        phase_shift=0.05,
        scheme="SPI",
        coherence_time=600.)
link_12.generate_entanglement(excitation_prob=0.03)
print(link_12.initial_fidelity_set)
chain = RelayChain(Link=link_12,
                   N=1)
chain.solve(swapping_quality=0.99)
print(chain.final_fidelity_without_decoherence)
print(chain.final_rate)

link_12 = HeraldedLink(
        detector_efficiency=1.,
        dark_count_prob=1.5e-6,
        transmission_without_fiber=0.08,
        distance=800./3,
        visibility=0.93,
        phase_shift=0.05,
        scheme="SPI",
        coherence_time=600.)
link_12.generate_entanglement(excitation_prob=0.03)
print(link_12.initial_fidelity_set)
chain = RelayChain(Link=link_12,
                   N=2)
chain.solve(swapping_quality=0.99)
print(chain.final_fidelity_without_decoherence)
print(chain.final_rate)

link_12 = HeraldedLink(
        detector_efficiency=1.,
        dark_count_prob=1.5e-6,
        transmission_without_fiber=0.08,
        distance=800./4,
        visibility=0.93,
        phase_shift=0.05,
        scheme="SPI",
        coherence_time=600.)
link_12.generate_entanglement(excitation_prob=0.03)
print(link_12.initial_fidelity_set)
chain = RelayChain(Link=link_12,
                   N=3)
chain.solve(swapping_quality=0.99)
print(chain.final_fidelity_without_decoherence)
print(chain.final_rate)

# 1000km
link_12 = HeraldedLink(
        detector_efficiency=1.,
        dark_count_prob=1.5e-6,
        transmission_without_fiber=0.08,
        distance=1000./2,
        visibility=0.93,
        phase_shift=0.05,
        scheme="SPI",
        coherence_time=600.)
link_12.generate_entanglement(excitation_prob=0.07)
print(link_12.initial_fidelity_set)
chain = RelayChain(Link=link_12,
                   N=1)
chain.solve(swapping_quality=0.99)
print(chain.final_fidelity_without_decoherence)
print(chain.final_rate)

link_12 = HeraldedLink(
        detector_efficiency=1.,
        dark_count_prob=1.5e-6,
        transmission_without_fiber=0.08,
        distance=1000./3,
        visibility=0.93,
        phase_shift=0.05,
        scheme="SPI",
        coherence_time=600.)
link_12.generate_entanglement(excitation_prob=0.07)
print(link_12.initial_fidelity_set)
chain = RelayChain(Link=link_12,
                   N=2)
chain.solve(swapping_quality=0.99)
print(chain.final_fidelity_without_decoherence)
print(chain.final_rate)

link_12 = HeraldedLink(
        detector_efficiency=1.,
        dark_count_prob=1.5e-6,
        transmission_without_fiber=0.08,
        distance=1000./4,
        visibility=0.93,
        phase_shift=0.05,
        scheme="SPI",
        coherence_time=600.)
link_12.generate_entanglement(excitation_prob=0.07)
print(link_12.initial_fidelity_set)
chain = RelayChain(Link=link_12,
                   N=3)
chain.solve(swapping_quality=0.99)
print(chain.final_fidelity_without_decoherence)
print(chain.final_rate)


link_12 = HeraldedLink(
        detector_efficiency=1.,
        dark_count_prob=1.5e-6,
        transmission_without_fiber=0.08,
        distance=1000./7,
        visibility=0.93,
        phase_shift=0.05,
        scheme="SPI",
        coherence_time=600.)
link_12.generate_entanglement(excitation_prob=0.07)
print(link_12.initial_fidelity_set)
chain = RelayChain(Link=link_12,
                   N=6)
chain.solve(swapping_quality=0.99)
print(chain.final_fidelity_without_decoherence)
print(chain.final_rate)









'''

gamma = 0.06
alpha = 0.03
alpha = alpha * (1-gamma) / (1 - alpha * gamma)

# 200km
link_12 = HeraldedLink(
        detector_efficiency=1.,
        dark_count_prob=1.5e-6,
        transmission_without_fiber=0.08,
        distance=200./2,
        visibility=0.93,
        phase_shift=0.05,
        scheme="SPI",
        coherence_time=600.)
link_12.generate_entanglement(excitation_prob=alpha)
print(link_12.initial_fidelity_set)
chain = RelayChain(Link=link_12,
                   N=1)
chain.solve(swapping_quality=0.99)
print(chain.final_fidelity_without_decoherence)
print(chain.final_rate)
'''