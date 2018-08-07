import nest
import numpy as np
import pdb

from brunel_alpha_nest import runBrunelNetwork as runBrunelNetworkAlpha

simtime = 1000
order = 500
eta = 2.0     # rel rate of external input
g = 5.0
dt = 0.1

NE = order*4
N = order*5

all_spikes = runBrunelNetworkAlpha(g=g, 
                              eta=eta, 
                              simtime=simtime, 
                              dt=dt,
                              order=order, 
                              simulator_name='NEST',
                              jnml_simulator=None)

# pdb.set_trace()