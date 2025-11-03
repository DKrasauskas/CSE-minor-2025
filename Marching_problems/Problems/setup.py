from Marching_problems.objects.function_definitions import *

# Consider the following problem:
#
#   ∂²u/∂t² − ∇ · (k ∇u) = f,  (x, y) ∈ Ω,  t ∈ (0, T]
#   u(x, y, t) = 0,             (x, y) ∈ ∂Ω
#   u(x, y, 0) = 0,             (x, y) ∈ Ω
#   ∂u/∂t (x, y, 0) = 0,       (x, y) ∈ Ω
#
# where the source term f is:
#   f(x, y, t) = sin(ω*t) * (
#       exp(-α*(x-x1)**2 - α*(y-y1)**2) +
#       exp(-α*(x-x2)**2 - α*(y-y2)**2) +
#       exp(-α*(x-x3)**2 - α*(y-y3)**2) +
#       exp(-α*(x-x4)**2 - α*(y-y4)**2)
#   )
#
# and the diffusion coefficient k(x, y) is piecewise:
#   k(x, y) = 0.1   if x < Lx/2 and y < Ly/2
#           = 0.4   if x < Lx/2 and y >= Ly/2
#           = 0.7   if x >= Lx/2 and y >= Ly/2
#           = 1.0   if x >= Lx/2 and y < Ly/2
#
# Parameters:
#   Ω = [0, Lx] × [0, Ly] with Lx = 10, Ly = 5
#   α = 40
#   (x1, y1) = (0.25*Lx, 0.25*Ly)
#   (x2, y2) = (0.25*Lx, 0.75*Ly)
#   (x3, y3) = (0.75*Lx, 0.75*Ly)
#   (x4, y4) = (0.75*Lx, 0.25*Ly)
#   ω = 4π


### Additionally, there is the diffusion problem:

# Consider the following problem:
#
#   ∂u/∂t − ∇ · (k ∇u) = f,  (x, y) ∈ Ω,  t ∈ (0, T]
#   u(x, y, t) = 0,           (x, y) ∈ ∂Ω
#   u(x, y, 0) = 0,           (x, y) ∈ Ω
#
# where the source term f is:
#   f(x, y) = exp(-α*(x-x1)**2 - α*(y-y1)**2)
#           + exp(-α*(x-x2)**2 - α*(y-y2)**2)
#           + exp(-α*(x-x3)**2 - α*(y-y3)**2)
#           + exp(-α*(x-x4)**2 - α*(y-y4)**2)
#
# and the diffusion coefficient k(x, y) is piecewise:
#   k(x, y) = 0.1   if x < Lx/2 and y < Ly/2
#           = 0.4   if x < Lx/2 and y >= Ly/2
#           = 0.7   if x >= Lx/2 and y >= Ly/2
#           = 1.0   if x >= Lx/2 and y < Ly/2
#
# Parameters:
#   Ω = [0, Lx] × [0, Ly] with Lx = 10, Ly = 5
#   α = 40
#   (x1, y1) = (0.25*Lx, 0.25*Ly)
#   (x2, y2) = (0.25*Lx, 0.75*Ly)
#   (x3, y3) = (0.75*Lx, 0.75*Ly)
#   (x4, y4) = (0.75*Lx, 0.25*Ly)

#both of these problems are implemented as Waves class and Diffusion class

L = (10, 5)
N = (200, 100)
Lx, Ly = L
Nx, Ny = N
p1 = np.array([ 0.25 * Lx, 0.25 * Ly])
p2 = np.array([ 0.25 * Lx, 0.75 * Ly])
p3 = np.array([0.75 * Lx, 0.75 * Ly])
p4 = np.array([ 0.75 *Lx, 0.25 * Ly])

launch_params_wave = (
    L,
    N,
    f_wave,
    (p1, p2, p3, p4),
    k,
    k_vectorized,
    None,
    0,
    4
)

launch_params_diffusion= (
    L,
    N,
    f,
    (p1, p2, p3, p4),
    k,
    k_vectorized,
    None,
)

