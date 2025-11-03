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
two substances — the ‘slow’ activator \( u \) and the ‘fast’ inhibitor \( v \) — 
leads to the emergence of regular patterns from noise. Such reactions occur, 
for instance, in animal skins and lead to the characteristic appearance of 
cheetahs and zebras.

---

## Model Equations

The process is described by the following system of non-linear coupled 
reaction-diffusion PDEs:

\[
\begin{aligned}
\frac{\partial u}{\partial t} &= D_u \, \Delta u + k \, (a - u + u^2 v), \quad &(1) \\
\frac{\partial v}{\partial t} &= D_v \, \Delta v + k \, (b - u^2 v), \quad &(2)
\end{aligned}
\]

for \((x, y) \in \Omega, \; t \in (0, T].\)

---

## Initial and Boundary Conditions

\[
\begin{aligned}
u(x, y, 0) &= u_0(x, y), \\
v(x, y, 0) &= v_0(x, y), \quad &(3)
\end{aligned}
\]

\[
\begin{aligned}
- D_u \, \nabla u \cdot n &= 0, \\
- D_v \, \nabla v \cdot n &= 0, \quad &(4)
\end{aligned}
\]

for \((x, y) \in \partial \Omega.\)

---

## Parameters

The rates of diffusion are determined by the corresponding diffusivity constants:

\[
D_u = 0.05, \quad D_v = 1.0
\]

The reaction constants are:

\[
k = 5, \quad a = 0.1305, \quad b = 0.7695
\]

---

## Initial Conditions

\[
\begin{aligned}
u_0(x, y) &= a + b + r(x, y), \\
v_0(x, y) &= \frac{b}{(a + b)^2}, \quad &(5)
\end{aligned}
\]

where \( r(x, y) \) is a small nonuniform perturbation in the concentration of 
the activator. All parameters are tuned to a regime where a pattern is expected 
to appear.

---

## Computational Domain

\[
\Omega = (0, 4) \times (0, 4)
\]

The pattern should be almost completely formed at:

\[
T = 20
\]



