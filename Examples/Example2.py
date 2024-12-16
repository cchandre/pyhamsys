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
potential = lambda x: -1 / np.sqrt(x**2 + 1)
psi0 = lambda x: np.exp(-x**2 / 50)
tf, tstep = 5e1, 1e-2
L, N = 100, 2**10

## Intermediate variables
xgrid = np.linspace(-L, L, N, endpoint=False, dtype=np.float64)
vgrid = potential(xgrid).astype(np.complex128)
dx = xgrid[1] - xgrid[0]

## Hamiltonian system and equations of motion
def y_dot(t, psi):
	psi_, _psi = np.roll(psi, 1), np.roll(psi, -1)
	lap_psi = (psi_ + _psi - 2 * psi) / dx**2
	return -1j * (-lap_psi / 2 + vgrid * psi) 

hs = HamSys(ndof=N, complex=True, y_dot=y_dot)

## Integration and plot of wavefunction
plt.ion()
fig, ax = plt.subplots()
h, = ax.plot(xgrid, np.abs(psi0(xgrid))**2, 'k', linewidth=2)
ax.set_title('$t = 0 $', loc='right', pad=20)
ax.set_xlabel('$x$')
ax.set_ylabel(r'$\vert\psi (x,t)\vert^2$')
ax.set_ylim((-0.1, 3.5))

def plot(t, psi):
	if (int(t / tstep) + 1) % 50 == 0:
		h.set_ydata(np.abs(psi)**2)
		ax.set_title(f'$t = {{{t:.1f}}}$', loc='right', pad=20)
		plt.pause(1e-4)
		
sol = solve_ivp_sympext(hs, (0, tf), psi0(xgrid), step=tstep, command=plot)

plt.ioff()
plt.show()