
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from renderer import *
from Objects.System import *

"""
The Schnakenberg model suggests a physical mechanism behind the emergence of 
Turing patterns in nature. The model shows that a chemical reaction between 
two substances — the ‘slow’ activator u and the ‘fast’ inhibitor v — leads to 
the emergence of regular patterns from noise. Such reactions occur, for instance, 
in animal skins and lead to the characteristic appearance of cheetahs and zebras.

The process is described by the following system of non-linear coupled 
reaction-diffusion PDEs:

    ∂u/∂t = Du * Δu + k * (a − u + u²v),        (1)
    ∂v/∂t = Dv * Δv + k * (b − u²v),            (2)

for (x, y) ∈ Ω, t ∈ (0, T];

with initial and boundary conditions:

    u(x, y, 0) = u₀(x, y),
    v(x, y, 0) = v₀(x, y),                      (3)

    −Du ∇u · n = 0,
    −Dv ∇v · n = 0,                            (4)
    for (x, y) ∈ ∂Ω.

The rates of diffusion are determined by the corresponding diffusivity constants:
    Du = 0.05
    Dv = 1.0

The reaction constants are:
    k = 5
    a = 0.1305
    b = 0.7695

The initial conditions are given by:
    u₀(x, y) = a + b + r(x, y),
    v₀(x, y) = b / (a + b)²,                   (5)

where r(x, y) is a small nonuniform perturbation in the concentration of 
the activator. All parameters are tuned to a regime where a pattern is expected 
to appear.

The computational domain is:
    Ω = (0, 4) × (0, 4)

The pattern should be almost completely formed at T = 20.
"""

INTEGRATION = {
    "EXPLICIT": 1,
    "IMPLICIT": 0
}


#initializing the list of parameters for the given PDE problem
params = (.1305, 0.7695, -0.05, -1.0, 4, 4, 100, 100)

#renderer for UI
ui = Renderer()

#creating the Coupled system
system = System(params, ui)

#if you want debugging info (like total wallclock etc), enable it by uncommenting the following line:
system.debug = True

#solve using either explicit or implicit integration model
system.solve(INTEGRATION["EXPLICIT"])
#uncomment for debug info
print(f"\033[33mwall clock time:\033[0m\n{system.wallclock}")
plt.show()

