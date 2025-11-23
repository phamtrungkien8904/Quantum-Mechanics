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


### Finite-Difference Derivation (GitHub‑Friendly Version)


To simulate 1D quantum systems numerically, we start from the time‑independent Schrödinger equation:


![eq](https://latex.codecogs.com/png.image?%5Cdpi%7B150%7D%20-%5Cfrac%7Bhbar%5E2%7D%7B2m%7D%20%5Cfrac%7Bd%5E2%20psi%7D%7Bdx%5E2%7D%20%2B%20V%28x%29%20psi%20%3D%20E%20\psi)


We discretize space into evenly spaced points x_j with spacing Δx. The wavefunction becomes values psi_j. The second derivative is approximated using the central finite‑difference formula:


![eq](https://latex.codecogs.com/png.image?%5Cdpi%7B150%7D%20%28psi_%7Bj%2B1%7D-2%20psi_j%20%2B%20psi_%7Bj-1%7D%29%20%2F%20%28%28Delta%20x%29%5E2%29)


Substituting into Schrödinger’s equation gives:


![eq](https://latex.codecogs.com/png.image?%5Cdpi%7B150%7D-%20%28hbar%5E2%2F%282m%29%29%20%28psi_%7Bj%2B1%7D-2%20psi_j%20%2B%20psi_%7Bj-1%7D%29%2F%28Delta%20x%29%5E2%20%2B%20V_j%20psi_j%20%3D%20E%20psi_j)


Define the constant λ:


![eq](https://latex.codecogs.com/png.image?%5Cdpi%7B150%7D%20lambda%20%3D%20hbar%5E2%20%2F%20%282m%20%28Delta%20x%29%5E2%29)


which yields the discrete equation:


![eq](https://latex.codecogs.com/png.image?%5Cdpi%7B150%7D%20-%20lambda%28psi_%7Bj%2B1%7D%20%2B%20psi_%7Bj-1%7D%29%20%2B%20%282lambda%20%2B%20V_j%29%20psi_j%20%3D%20E%20psi_j)


Collecting all values psi_j into a vector gives the matrix eigenvalue problem:


![eq](https://latex.codecogs.com/png.image?%5Cdpi%7B150%7D%20H%20Psi%20%3D%20E%20Psi)


The Hamiltonian matrix H is tridiagonal:


![matrix](https://latex.codecogs.com/png.image?%5Cdpi%7B150%7D%20%5Cbegin%7Bpmatrix%7D2lambda%20%2B%20V_1%20%26%20-lambda%20%26%200%20%26%20%5Ccdots%20%5C%5C-lambda%20%26%202lambda%20%2B%20V_2%20%26%20-lambda%20%26%20%5Ccdots%20%5C%5C0%20%26%20-lambda%20%26%202lambda%20%2B%20V_3%20%26%20%5Ccdots%20%5C%5C%5Cvdots%20%26%20%5Cvdots%20%26%20%5Cvdots%20%26%20%5Cddots%5Cend%7Bpmatrix%7D)


Boundary conditions for the infinite square well require:


![eq](https://latex.codecogs.com/png.image?%5Cdpi%7B150%7D%20psi_0%20%3D%20psi_%7BN%2B1%7D%20%3D%200)


Thus, only psi_1 through psi_N appear in the matrix. Solving the eigenvalue problem yields:
- energy levels E_n
- stationary wavefunctions psi_n
