#
# BSD 2-Clause License
#
# Copyright (c) 2024, Cristel Chandre
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from pyhamsys import HamSys

## Parameters
epsilon = 0.027					# parameter of the Hamiltonian system
timestep = 2 * np.pi / 50		# integration timestep 
N = 50  						# number of trajectories
nf = 1000						# number of points on the Poincaré section per trajectory

## Hamiltonian system
def hamiltonian(q, p, t):
	return sum([p_**2 / 2 + epsilon * (sp.cos(q_ - t) + sp.cos(q_)) for q_, p_ in zip(q, p)])

hs = HamSys(ndof=N)
hs.compute_vector_field(hamiltonian)
	
## Initial conditions
x0 = 2 * np.pi * np.random.random(N)
p0 = np.random.random(N)
y0 = np.concatenate((x0, p0), axis=None)

## Integration
sol = hs.integrate(y0, 2 * np.pi * np.arange(nf + 1), extension=True, timestep=timestep)

## Plot of the Poincaré section
fig, ax = plt.subplots()
x, p = np.split(sol.y, 2)
ax.plot(x % (2 * np.pi), p, '.b')
ax.set_xlim((0, 2 * np.pi))
plt.show()