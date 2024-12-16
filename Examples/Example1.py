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
import matplotlib.pyplot as plt
from pyhamsys import HamSys, solve_ivp_sympext

## Parameters
epsilon = 0.027				# parameter of the Hamiltonian system
tstep = 2 * np.pi / 10		# integration timestep 
N = 200  					# number of trajectories
nf = 1500 					# number of points on the Poincaré section per trajectory

## Hamiltonian systems and equations of motion 
def y_dot(t, y):
	y_ = np.split(y, 2)
	return np.concatenate((y_[1], epsilon * (np.sin(y_[0] - t) + np.sin(y_[0]))), axis=None)

def hamiltonian(t, y):
	y_ = np.split(y, 2)
	return np.sum(y_[1]**2 / 2 + epsilon * (np.cos(y_[0] - t) + np.cos(y_[0])))

def k_dot(t, y):
	y_ = np.split(y, 2)
	return -epsilon * np.sum(np.sin(y_[0] - t))

hs = HamSys(ndof=N + 0.5, y_dot=y_dot, k_dot=k_dot, hamiltonian=hamiltonian)
	
x0 = 2 * np.pi * np.random.random(N)
p0 = np.random.random(N)
y0 = np.concatenate((x0, p0), axis=None)

## Integration
sol = solve_ivp_sympext(hs, (0, 2 * np.pi * nf), y0, step=tstep, t_eval=2 * np.pi * np.arange(nf + 1), check_energy=True)
print(sol.err)

## Plot of the Poincaré section
fig, ax = plt.subplots()
x, p = np.split(sol.y, 2)
ax.plot(x % (2 * np.pi), p, '.b')
ax.set_xlim((0, 2 * np.pi))
plt.show()