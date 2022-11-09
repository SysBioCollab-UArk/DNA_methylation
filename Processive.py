from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

Model()
'''
Molecules

H,H, hemimethylated CpG, is hemethylated sites( sites that have not been remethylated)the empty pentagons
M, means methylated sites, represented as red pentagons
SAM molecule, represented as S, S-adeno- sylmethionine (SAM)
SAH,unmethylated form of SAM, represented by Q IN THE PICTURE, S-adenosyl homocysteine (SAH)


DISTRIBUTIVE MODEL
(STEP 2):DNMT1 binds to a hemimethylated CpG site,
(STEP 2): incorporates SAM
(STEP 3) catalizes methylation
Methylation release both DNMT1 and SAH occur in a single step

PROCESSIVE MODEL
After DNMT1 can diffuse towards  its immediate neighbor sites either upstream(U)
or downstream(D)
Assumes that DNMT1 can remain bound to DNA after performing the catalytic step 
and reach the subsequent hemimethylated CpG sites by diffusing along DNA

The first two steps are identical to Distributive Model
However, Enzyme remains bound to DNA molecule onto the recenlty methylated site 
forming the EM complex


COLLABORATIVE MODEL
Once DNMT1 is bound to a hemimethylated CpG,a second DNMT1 molecule can be recruited 
onto nearby CpGsites

'''
# Monomers
Monomer('DNMT1', ['h', 'sam'])
Monomer('SAM', ['dnmt1'])
Monomer('H', ['dnmt1'])
Monomer('M',['dnmt1'])
Monomer('SAH')

# Initials Conditions
Parameter('DNMT1_0', 100)
Parameter('H_0', 200)
Parameter('SAM_0', 150)
Initial(DNMT1(h=None, sam=None), DNMT1_0)
Initial(H(dnmt1=None), H_0)
Initial(SAM(dnmt1=None), SAM_0)

# Parameter
Parameter('kf1', 100)
Parameter('kr1', 1)
Parameter('kf2', 100 )
Parameter('kr2', 100 )
Parameter('k3', 1)
Parameter('kDif_U', 100)
Parameter('kDif_D', 100)
Parameter('koff', 100)

# Rules

Rule('DNMT1_binds_H', DNMT1(h=None, sam=None) + H(dnmt1=None) | DNMT1(h=1, sam=None) % H(dnmt1=1), kf1, kr1)
Rule('SAM_binds_DNMT1_H', SAM(dnmt1=None) + DNMT1(h=1, sam= None) % H(dnmt1=1) | SAM(dnmt1=2) % DNMT1(h=1, sam=2) % H(dnmt1=1), kf2, kr2)
Rule('DNMT1_is_methylated', DNMT1(h=1, sam=2) % SAM(dnmt1=2) % H(dnmt1=1) >> DNMT1(h=1,sam=None) % M(dnmt1=1) + SAH() , k3)
Rule('DNMT1_diffuses_Upstream', DNMT1(h=1, sam=None) % M(dnmt1=1) + H(dnmt1=None) >> DNMT1(h=1, sam=None) % H(dnmt1=1) + M(dnmt1=None), kDif_U)
Rule('DNMT1_diffuses_downstream', DNMT1(h=1, sam=None) % M(dnmt1=1) + H(dnmt1=None) >> DNMT1(h=1, sam=None) % H(dnmt1=1) + M(dnmt1=None), kDif_D)
Rule('DNMT1_Unbinds_M', DNMT1(h=1, sam=None) % M(dnmt1=1) >> DNMT1(h=None, sam=None)+ M(dnmt1=None), koff)

# Observables
Observable('H_total', H())
Observable('H_bound', H(dnmt1= ANY))
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
