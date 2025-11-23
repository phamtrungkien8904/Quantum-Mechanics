# Quantum-Mechanics
A small and generalized Python repository that simulates several 1D quantum system using the finite-difference method (FDM) and linear algebra to compute eigenvalues and eigenvectors of the time-independent Schrödinger equation.

# Quantum-Mechanics


Simulate 1D quantum systems using the finite-difference method and linear algebra to solve for eigenvalues (energies) and eigenvectors (stationary states).


## Included scripts


- `infinite-well.py` — particle in an infinite potential well (box). Computes eigenstates and makes an animation of the first few stationary states.
- `tunneling.py` — wavefunction stationary states for a finite barrier and an animation showing a wavepacket encountering a barrier (tunneling/transmission/reflection).
- `oscilator.py` — quantum harmonic oscillator eigenstates and animation of a superposition.
- `free-fall.py` — particle in a linear gravitational-like potential (V(x)=m g x) and eigenstates/animation.


## Method


### Finite-Difference Derivation (Full Explanation)


To simulate 1D quantum systems numerically, we start from the time-independent Schrödinger equation:


![equation](https://latex.codecogs.com/svg.image?-\frac{\hbar^2}{2m}\frac{d^2\psi}{dx^2}&plus;V(x)\,\psi(x)=E\,\psi(x).)

Space is discretized into evenly spaced points x_j with spacing Δx. The wavefunction becomes values ψ_j = ψ(x_j). The second derivative is approximated using the central finite-difference formula:

```math
\frac{d^2\psi}{dx^2} \approx \frac{\psi_{j+1} - 2\psi_j + \psi_{j-1}}{(\Delta x)^2}. \]
```

Substituting into the Schrödinger equation yields:


\[ -\frac{\hbar^2}{2m} \frac{\psi_{j+1} - 2\psi_j + \psi_{j-1}}{(\Delta x)^2} + V_j \psi_j = E \psi_j. \]


Define:


\[ \lambda = \frac{\hbar^2}{2m(\Delta x)^2}. \]


Then the discrete equation becomes:


\[ -\lambda (\psi_{j+1} + \psi_{j-1}) + (2\lambda + V_j)\psi_j = E\psi_j. \]


This forms a linear algebra eigenvalue problem:


\[ H \Psi = E \Psi. \]


The Hamiltonian matrix H is tridiagonal:


\[
H = \begin{pmatrix}
2\lambda + V_1 & -\lambda & 0 & \cdots \\
-\lambda & 2\lambda + V_2 & -\lambda & \cdots \\
0 & -\lambda & 2\lambda + V_3 & \cdots \\
\vdots & \vdots & \vdots & \ddots
\end{pmatrix}. \]


Boundary conditions determine the top and bottom rows. For the infinite square well, ψ_0 = ψ_{N+1} = 0, so only ψ_1 ... ψ_N are included. Solving the eigenvalue problem yields numerical approximations of the energy levels E_n and stationary states ψ_n.


This project uses the standard finite-difference discretization of the time-independent Schr\u00f6dinger equation. Replacing derivatives with central differences converts


\[ -\frac{\hbar^2}{2m} \frac{d^2\psi}{dx^2} + V(x)\,\psi = E\,\psi \]


into the discrete form


\[
-\lambda(\psi_{j+1} + \psi_{j-1}) + (2\lambda + V_j)\psi_j = E\psi_j,
\]


where \lambda = \frac{\hbar^2}{2m(\Delta x)^2}. This yields the matrix eigenvalue equation


\[ H\Psi = E\Psi \]


with a tridiagonal Hamiltonian. Eigenvalues approximate allowed energies; eigenvectors approximate stationary wavefunctions.


All scripts use the finite-difference discretization of the 1D time-independent Schrödinger equation


\[-\frac{\hbar^2}{2m}\frac{d^2}{dx^2} + V(x)\] \psi(x) = E \psi(x)


The second derivative is approximated with central finite differences on a uniform grid, leading to a symmetric tridiagonal Hamiltonian matrix. We solve for eigenvalues and eigenvectors using `scipy.linalg.eigh`.


## Requirements


Python 3.9+ and the packages listed in `requirements.txt`.


## Install


```bash
python -m venv venv
source venv/bin/activate # or `venv\\Scripts\\activate` on Windows
pip install -r requirements.txt
