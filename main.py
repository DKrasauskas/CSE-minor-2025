from Problems.setup import *
from objects.function_definitions import *
from Problems.waves.waves import *

w = Waves(L, N, launch_params_wave)
d = Diffusion(L, N, launch_params_diffusion)
w.get()

