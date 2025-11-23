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


### Finite-Difference Derivation (GitHub-Friendly Version)


To simulate 1D quantum systems numerically, we start from the time-independent Schrödinger equation:


![eq](https://latex.codecogs.com/png.image?%5Cdpi%7B150%7D&fg=white%20-%5Cfrac%7Bhbar%5E2%7D%7B2m%7D%20%5Cfrac%7Bd%5E2%5Cpsi%28x%29%7D%7Bdx%5E2%7D%20%2B%20V%28x%29%5Cpsi%28x%29%20%3D%20E%5Cpsi%28x%29)


We discretize space into evenly spaced points x_j with spacing Δx. The wavefunction becomes values ψ(x_j). The second derivative is approximated using the central finite-difference formula:


![eq](https://latex.codecogs.com/png.image?%5Cdpi%7B150%7D&fg=white%20%28%5Cpsi%28x_%7Bj%2B1%7D%29-2%5Cpsi%28x_j%29%2B%5Cpsi%28x_%7Bj-1%7D%29%29%2F%28%28%5CDelta%20x%29%5E2%29)


Substituting into Schrödinger’s equation gives:


![eq](https://latex.codecogs.com/png.image?%5Cdpi%7B150%7D&fg=white%20-%20%28hbar%5E2%2F%282m%29%29%20%28%5Cpsi%28x_%7Bj%2B1%7D%29-2%5Cpsi%28x_j%29%2B%5Cpsi%28x_%7Bj-1%7D%29%29%2F%28%5CDelta%20x%29%5E2%20%2B%20V_j%5Cpsi%28x_j%29%20%3D%20E%5Cpsi%28x_j%29)


Define the constant λ:


![eq](https://latex.codecogs.com/png.image?%5Cdpi%7B150%7D&fg=white%20%5Clambda%20%3D%20hbar%5E2%2F%282m%28%5CDelta%20x%29%5E2%29)


which yields the discrete equation:


![eq](https://latex.codecogs.com/png.image?%5Cdpi%7B150%7D&fg=white%20-%20%5Clambda%28%5Cpsi%28x_%7Bj%2B1%7D%29%20%2B%20%5Cpsi%28x_%7Bj-1%7D%29%29%20%2B%20%282%5Clambda%20%2B%20V_j%29%5Cpsi%28x_j%29%20%3D%20E%5Cpsi%28x_j%29)


The matrix eigenvalue form is:


![eq](https://latex.codecogs.com/png.image?%5Cdpi%7B150%7D&fg=white%20H%20%5CPsi%28x%29%20%3D%20E%20%5CPsi%28x%29)


The Hamiltonian matrix H is tridiagonal:


![matrix](https://latex.codecogs.com/png.image?%5Cdpi%7B150%7D&fg=white%20%5Cbegin%7Bpmatrix%7D2%5Clambda%2BV_1%20%26-lambda%20%260%20%26%5Ccdots%5C%5C-lambda%20%262lambda%2BV_2%20%26-lambda%20%26%5Ccdots%5C%5C0%20%26-lambda%20%262lambda%2BV_3%20%26%5Ccdots%5C%5C%5Cvdots%20%26%5Cvdots%20%26%5Cvdots%20%26%5Cddots%5Cend%7Bpmatrix%7D)


Boundary conditions for the infinite well:


![eq](https://latex.codecogs.com/png.image?%5Cdpi%7B150%7D&fg=white%20%5Cpsi%280%29%3D%5Cpsi%28L%29%3D0)


Solving yields numerical approximations to energy levels E_n and wavefunctions ψ_n(x).
