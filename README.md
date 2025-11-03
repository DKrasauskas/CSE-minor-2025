# *Numerical methods for PDE's*


This repo contains all of the code I used in my CSE minor course  Numerical methods for PDEs. It contains a set of Finite-Volume and Finitie-Difference problems. 


<table>
  <tr>
     <td>
      <img src="https://github.com/user-attachments/assets/3d742229-25ea-42c3-8f94-5d342aace95d" width="100%" >
        <figcaption align="center"> Schnakenberg at T = 0 </figcaption>
    </td>
    <td>
      <img src="https://github.com/user-attachments/assets/1ddd1714-43da-46f4-a1d7-8038abc5e7f7" width="100%">
       <figcaption align="center"> Schnakenberg at T = 20 </figcaption>
    </td>
  </tr>
</table>

# The Schnakenberg Model

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

# Boundary-Value Problem

Consider the following boundary value problem:

    −∇ · (k∇u) = f,          (x, y) ∈ Ω = (0, 10) × (0, 5)
    u(x, y) = 0,              (x, y) ∈ ∂Ω

where

    f(x, y) = 
        Σ (i = 1 to 9) Σ (j = 1 to 4)
        e^(−α(x−i)²−α(y−j)²),

with

    α = 40,   (x, y) ∈ Ω.                          (1)

Here, Ω = [0, 10] × [0, 5] is the rectangle with the corners (0, 0), (10, 0), (10, 5), and (0, 5).

---

## Finite-Difference Method

For this method we shall only consider the homogeneous coefficient function:

    k(x, y) = 1,     (x, y) ∈ Ω.

---

## Finite-Volume Method

For this method we shall consider both the homogeneous coefficient function (2) and the following inhomogeneous coefficient function:

    k(x, y) = 1 + 0.1(x + y + xy),     (x, y) ∈ Ω.


The following results are obtained:
<table>
  <tr>
     <td>
      <img src="https://github.com/user-attachments/assets/8d96b60d-fe2d-42c4-acd4-cef5caf39f21"  width="100%" >
        <figcaption align="center"> Schnakenberg at T = 0 </figcaption>
    </td>
    <td>
      <img  src="https://github.com/user-attachments/assets/682c569d-398b-4e64-909c-20dac8f1b50a" width="100%">
       <figcaption align="center"> Schnakenberg at T = 20 </figcaption>
    </td>
  </tr>
</table>

The following results are obtained:
<table>
  <tr>
     <td>
      <img src="https://github.com/user-attachments/assets/9dc67b2b-6564-4166-a8dc-74f45e5f1d9b"  width="100%" >
        <figcaption align="center"> Schnakenberg at T = 0 </figcaption>
    </td>
    <td>
      <img src="https://github.com/user-attachments/assets/f3cf0b4f-1758-4a9f-a5f4-2ccaced42184" width="100%">
       <figcaption align="center"> Schnakenberg at T = 20 </figcaption>
    </td>
  </tr>
</table>

The pattern should be almost completely formed at T = 20.



