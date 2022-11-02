from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

Model()

# Monomers
Monomer('DNMT1', ['h', 'sam'])
Monomer('H', ['dnmt1'])
Monomer('SAM', ['dnmt1'])
Monomer('M')
Monomer('SAH')

# Initial Conditions
Parameter('DNMT1_0', 10)
Parameter('H_0', 2000)
Parameter('SAM_0', 1000)
Initial(DNMT1(h=None, sam=None), DNMT1_0)
Initial(H(dnmt1=None), H_0)
Initial(SAM(dnmt1=None), SAM_0)

# Rules
Parameter('kf1', 100)
Parameter('kr1', 1)
Parameter('kf2', 100)
Parameter('kr2', 1)
Parameter('k3', 100)
Rule('DNMT1_binds_H', DNMT1(h=None, sam=None) + H(dnmt1=None) | DNMT1(h=1, sam=None) % H(dnmt1=1), kf1, kr1)

Rule('SAM_binds_DNMT1_H', SAM(dnmt1=None) + DNMT1(h=1, sam=None) % H(dnmt1=1) |
     SAM(dnmt1=2) % DNMT1(h=1, sam=2) % H(dnmt1=1), kf2, kr2)

Rule('DNA_methylation', SAM(dnmt1=2) % DNMT1(h=1, sam=2) % H(dnmt1=1) >> DNMT1(h=None, sam=None) + SAH() + M(), k3)

# Observables
Observable('H_tot', H())
Observable('H_bound', H(dnmt1=ANY))
Observable('M_tot', M())
Observable('DNMT1_unbound', DNMT1(h=None, sam=None))
Observable('SAM_tot', SAM())

# Actions
tspan = np.linspace(0, 2, 201)
sim = ScipyOdeSimulator(model, tspan, verbose=True)
output = sim.run()

for obs in model.observables:
    plt.plot(tspan, output.observables[obs.name], lw=2, label=obs.name)
plt.legend(loc=0)
plt.xlabel('time')
plt.ylabel('amount')

plt.tight_layout()
plt.show()
