from pysb import*
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

#Collaborative Model
Model()

# Monomers
Monomer('DNMT1', ['h', 'sam'])
Monomer('SAM', ['dnmt1'])
Monomer('H', ['dnmt1'])
Monomer('M')
Monomer('SAH')

# Initial Conditions
Parameter('DNMT1_0', 100)
Parameter('DNMT1_A', 100)
Parameter('H_0', 200)
Parameter('SAM_0', 150)
Initial(DNMT1(h=None,sam=None), DNMT1_0)
Initial(H(dnmt1=None), H_0)
Initial(SAM(dnmt1=None), SAM_0)

# Rules
Parameter('kf1', 100)
Parameter('kr1', 0)
Parameter('kf2', 100)
Parameter('kr2', 0)
Parameter('k3', 100)
Parameter('krec_U', 100)
Parameter('krec_D', 100)

#Rule
Rule('DNMT1_binds_H', DNMT1(h=None, sam=None) + H(dnmt1=None) | DNMT1(h=1, sam=None) % H(dnmt1=1), kf1, kr1)

Rule('SAM_binds_DNMT1_H', SAM(dnmt1=None) + DNMT1(h=1, sam=None) % H(dnmt1=1) | SAM(dnmt1=2) % DNMT1(h=1, sam=2) % H(dnmt1=1), kf2, kr2)

Rule('DNA_methylation', SAM(dnmt1=2) % DNMT1(h=1, sam=2) % H(dnmt1=1) >> DNMT1(h=None, sam=None) + SAH() + M(), k3)

'''
here sites are interdependent through a phenomenological mechanism of “collaboration” between 
enzyme molecules: after DNMT1 binds, a second enzyme can be recruited to any nearby CpG site 
upstream or downstream the original site (not necessarily the contiguous). The stochastic 
pro- pensity for each recruitment reaction kRec, notwithstanding, is indeed a function of the
distance of a neighboring hemimethylated site to the recruiting site according to:

'''
#after DNMT1 binds, a second enzyme can be recruited to any nearby CpG site upstream or
# downstream the original site (not necessarily the contiguous).

Rule('DNMT1_krec_U', DNMT1(h=1, sam=None) % H(dnmt1=1) + DNMT1(h=None, sam=None) + H(dnmt1=None) >>
     DNMT1(h=1, sam=None) % H(dnmt1=1) + DNMT1(h=2, sam=None) % H(dnmt1=2), krec_U)

Rule('DNMT1_krec_D', DNMT1(h=1, sam=None) % H(dnmt1=1) + DNMT1(h=None, sam=None) + H(dnmt1=None) >>
     DNMT1(h=1, sam=None) % H(dnmt1=1) + DNMT1(h=2, sam=None) % H(dnmt1=2), krec_D)


#Observables
Observable('H_total', H())
Observable('H_bound', H(dnmt1=ANY))
Observable('M_total', M())
Observable('DNMT1_unbound', DNMT1(h=None, sam=None))
Observable('SAM_total', SAM())

#Actions
tspan = np.linspace(0, 2, 100)
sim =ScipyOdeSimulator(model, tspan, verbose=True)
output = sim.run()

for obs in model.observables:
    plt.plot(tspan, output.observables[obs.name], lw=2, label=obs.name)
    plt.legend(loc=0)
    plt.xlabel('time')
    plt.ylabel('amount')

    plt.tight_layout()
    plt.show()
