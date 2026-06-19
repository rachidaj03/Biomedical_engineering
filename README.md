# Discrétisation 1 — Numerical Methods for Biomedical Engineering (S5)

Project on **finite-difference discretization of PDEs**, applied to modeling oxygen (O₂) diffusion/transport in biological tissue (alveolar–capillary exchange). The work progresses from simple 1D steady-state problems to 2D steady-state and transient problems on regular and arbitrary (irregular) domains.

## Structure

### `Ordre 1/` — 1D stationary, first-order scheme
`1D_stationnaire_ordre_1.ipynb`: solves `u'(x) + λu(x) = 0` on `[0, 1]` (λ = -3, u(0) = 9) with a first-order finite-difference scheme (tridiagonal system, n = 100 points). Compares the numerical solution to the exact solution `u(x) = 9·exp(-3x)`.

### `Ordre 2/` — 1D stationary, second-order scheme
`1D_stationnaire_ordre_2.ipynb`: solves `u''(x) + λu(x) = 0` (λ = 90) with a second-order central-difference scheme (n = 1000 points), and compares against the exact solution `u(x) = α·exp(√λ·x) + β·exp(-√λ·x)`.

### `1D temporelle/` — 1D transient schemes
Four notebooks comparing time-stepping schemes for advection (`∂u/∂t + c·∂u/∂x = 0`) and diffusion (`∂u/∂t = c·∂²u/∂x²`), initial condition `u(x,0) = cos(x)`:
- `Schema_2_diffusion.ipynb` — implicit (backward) Euler
- `Schema_2_advection.ipynb` — explicit upwind scheme
- `Schema_3_diffusion.ipynb` — Crank-Nicolson (2nd-order in time)
- `Schema_3_advection.ipynb` — Lax-Wendroff-style 3-level scheme (leapfrog)

Used to study stability (CFL/Courant number) and accuracy trade-offs between explicit and implicit methods.

### `2D_non_temporelle/` — 2D steady-state O₂ transport
2D Poisson/Helmholtz-type problem `∇²u + λu = 0` on a rectangular domain, with Dirichlet (top), Neumann (sides), and Robin (bottom, absorption) boundary conditions.
- `Final.py` — base rectangular-domain solver (generic parameters)
- `Final_2.py` — same model with real physiological parameters (O₂ diffusivity, membrane permeability, alveolar/venous partial pressures), including inlet/outlet flux balance
- `dom_quelc.py`, `domaine-2.py` — extend the solver to irregular geometry (bronchi channels, alveolar absorption zones) via a mesh activation matrix encoding interior/Dirichlet/Neumann/Robin/corner points
- `absorption_imparfaite.py`, `analytique.py`, `frequentiel.py`, `temporel.py`, `temporelle_final.py`, `P1.py`, `fichier.py`, `tring.py`, `essai.py` — exploratory/intermediate variants (imperfect absorption BC, analytical comparison, frequency-domain analysis, transient extensions)

### `Domaine quleconque/` — arbitrary domain handling
Earlier/parallel exploration of solving on non-rectangular domains: `domaine.py`, `domaine-2.py`, `code.py`, `test.py`, `test_2.py`, `Domaine.ipynb`. `output.csv` holds an exported solution grid.

### Root files
- `Test.ipynb` — 2D transient diffusion (implicit Euler) with interactive 3D Tkinter visualization
- `Or.py` — unrelated scratch file (DFS maze generator), not part of the PDE project
- `myAnimation.gif` — animation of a simulation's time evolution

## Key concepts covered
- Finite-difference discretization (1st & 2nd order) of ODEs/PDEs
- Dirichlet, Neumann, and Robin boundary conditions
- Explicit vs. implicit time-stepping, stability (CFL condition), Crank-Nicolson
- 2D Laplacian/Helmholtz assembly via grid-to-vector index mapping
- Irregular-domain handling via mesh activation/masking
- Application: O₂ diffusion through tissue with membrane permeability and capillary absorption

## Stack
Python (Jupyter notebooks + scripts), NumPy, SciPy (`linalg`), Matplotlib, Tkinter (3D interactive plots).
