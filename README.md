# Quantum-Mechanics
A small and generalized Python repository that simulates several 1D quantum systems with arbitrary potential distribution and initial wave function using the finite-difference method (FDM) and linear algebra to compute eigenvalues and eigenvectors of the time-independent Schrödinger equation.




## Included scripts


- `infinite-well.py` — particle in an infinite potential well (box). Computes eigenstates and makes an animation of the first few stationary states.
- `tunneling.py` — wavefunction stationary states for a finite barrier with potential $`V_01$ and an animation showing a wavepacket encountering a barrier.
- `oscilator.py` — quantum harmonic oscillator eigenstates and animation of a superposition.
- `free-fall.py` — particle in a linear gravitational-like potential ($`V(x)=m g x`$) and eigenstates/animation.


## Method


### Finite-Difference Derivation (Full Explanation)


To simulate 1D quantum systems numerically, we start from the time-independent Schrödinger equation:


```math
-\frac{\hbar^2}{2m} \frac{d^2\psi}{dx^2} + V(x)\,\psi(x) = E\,\psi(x).
```

Space is discretized into evenly spaced points $`x_j`$ with spacing $`\Delta`$. The wavefunction becomes values $`\psi_j = \psi(x_j)`$. The second derivative is approximated using the central finite-difference formula:

```math
\frac{d^2\psi}{dx^2} \approx \frac{\psi_{j+1} - 2\psi_j + \psi_{j-1}}{(\Delta x)^2}. 
```

Substituting into the Schrödinger equation yields:


```math
-\frac{\hbar^2}{2m} \frac{\psi_{j+1} - 2\psi_j + \psi_{j-1}}{(\Delta x)^2} + V_j \psi_j = E \psi_j. 
```

Define:

```math
\lambda = \frac{\hbar^2}{2m(\Delta x)^2}.
```

Then the discrete equation becomes:

```math
-\lambda (\psi_{j+1} + \psi_{j-1}) + (2\lambda + V_j)\psi_j = E\psi_j.
```

This forms a linear algebra eigenvalue problem:

```math
\hat{H} \psi_n = E_n \psi_n.
```

The Hamiltonian matrix $`H`$ is tridiagonal:


```math
\hat{H} = \begin{pmatrix}
2\lambda + V_1 & -\lambda & 0 & \cdots \\
-\lambda & 2\lambda + V_2 & -\lambda & \cdots \\
0 & -\lambda & 2\lambda + V_3 & \cdots \\
\vdots & \vdots & \vdots & \ddots
\end{pmatrix}. 
```

Boundary conditions determine the top and bottom rows. For the infinite square well, $`\psi_0 = \psi_{N+1}`$ = 0, so only $`\psi_1 \ldots \psi_N`$ are included. Solving the eigenvalue problem yields numerical approximations of the energy levels $`E_n`$ and stationary states $`\psi_n`$.


The second derivative is approximated with central finite differences on a uniform grid, leading to a symmetric tridiagonal Hamiltonian matrix. We solve for eigenvalues and eigenvectors using `np.linalg.eigh`.

After that we can find the general solution of the wave function $Psi(x,t)$. Firstly, we find the expansion coefficient for corresponding eigen wavefunction: c_n = \int psi_n^*(x) \Psi(x,0) dx. Finally, \Psi(x,t) = \sum c_n \psi_n(x) e^{-i E_n t/\hbar}.
