# Fourier–Collocation Framework for QPO Family Construction in the CR3BP

**Authors:** Evan Lantier, Roshan Eapen
**Affiliation:** Department of Aerospace Engineering, The Pennsylvania State University




---

## Overview

This repository contains the MATLAB implementation of a Fourier–collocation framework for constructing families of quasi-periodic orbits (QPOs) in the Earth–Moon Circular Restricted Three-Body Problem (CR3BP).

Starting from a single reference periodic orbit, the framework:

1. Differentially corrects the periodic orbit to machine precision
2. Extracts fundamental frequencies via Floquet theory and the NAFF algorithm
3. Constructs an initial Fourier-based torus approximation seeded from the periodic orbit
4. Refines the torus to enforce dynamical consistency with the true CR3BP flow
5. Continues the orbit into a full QPO family via arclength continuation in Fourier coefficient space
6. Systematically identifies and adds missing harmonic combinations to improve accuracy

The result is a semi-analytical torus representation — expressed directly in terms of fundamental frequencies and Fourier coefficients — that can be evaluated at any point, continued into families, and quantitatively validated through phase-point drift analysis.

---

## Background

The CR3BP is a Hamiltonian system governing the motion of a massless particle under the gravitational influence of two primaries (Earth and Moon) in a rotating reference frame. Near the collinear Lagrange points (L1, L2), the linearized dynamics admit a saddle × center × center structure, giving rise to families of periodic and quasi-periodic orbits organized around two independent frequencies.

A quasi-periodic orbit on a 2-torus is parameterized by two angle variables **Q** = [Q₁, Q₂] evolving at rates **Ω** = [ω₀, ωₚ], and its Cartesian coordinates are represented as a truncated Fourier series:

```
q_i(Q) = Σ_j [ A_j · cos(j · Q) + B_j · sin(j · Q) ]
```

where **j** ∈ ℤ² are integer combination vectors. This framework constructs and refines the coefficient arrays {A, B} and the frequency vector **Ω** to satisfy the CR3BP equations of motion.

---

## Repository Structure

```
.
├── Fourier_Collocation_Template.m   # Main script + all functions (self-contained)
└── README.md
```

The implementation is contained in a single self-documented MATLAB file. All functions are defined inline at the bottom of the script.

### Script Sections

| Section | Description |
|---|---|
| 1 | Periodic orbit initialization and differential correction |
| 2 | Monodromy matrix, Floquet multipliers, Poincaré exponents |
| 3 | PO propagation and plotting |
| 4 | Spectral decomposition of the PO via least-squares Fourier fit |
| 5 | QPO frequency grid construction (`Freq_Comb`) |
| 6 | Initial QPO refinement at ε = ds from the reference PO |
| 7 | QPO family continuation (arclength in Fourier space) |
| 8 | Systematic harmonic addition |
| 9 | Trajectory plots (second-order and added-harmonic families) |
| 10 | Phase-point drift validation and violin plots |

### Functions

| Function | Purpose |
|---|---|
| `CR3BP` | CR3BP equations of motion |
| `CR3BP_STM` | CR3BP EOM + State Transition Matrix (42-state) |
| `JacobiConstant` | Compute Jacobi constant C for a state array |
| `Freq_Comb` | Generate integer frequency combination grid |
| `QPO_fun_first` | Objective: Jacobi constant residual (initial QPO) |
| `QPO_con_first` | Constraints: collocation + center manifold displacement |
| `QPO_fun` | Objective: Jacobi constant residual (continuation) |
| `QPO_con` | Constraints: collocation + arclength + phase conditions |
| `QPO_con_added` | Constraints: collocation + IC pinning (harmonic addition) |
| `QPO_Add_Harmonics` | Integer relation fit via CVX/Gurobi to find missing harmonics |
| `NAFF` | Numerical Analysis of Fundamental Frequencies algorithm |
| `Q_to_q` / `Q_to_v` | Evaluate torus position/velocity from angle coordinates |
| `ConstructOrbit` | Evaluate Fourier orbit from single-frequency basis |
| `ExtractFreq` | FFT peak extraction for NAFF initialization |
| `ShootingParams` | Inner products for NAFF Newton iteration |
| `JacobianComponents` | 3×3 Jacobian for NAFF frequency refinement |
| `ReduceSignal` / `ReconSig` | NAFF signal manipulation utilities |

---

## Requirements

### MATLAB
- MATLAB R2020b or later (tested)
- Optimization Toolbox — `lsqnonlin` with `interior-point` algorithm
- Parallel Computing Toolbox — recommended for `UseParallel: true` in `lsqnonlin`

### Third-Party (Required)
- **CVX** (v2.1+) — convex optimization modeling framework
  - https://cvxr.com/cvx/download/
- **Gurobi** (v9.0+) — mixed-integer solver called by CVX
  - https://www.gurobi.com (free academic license available)

### Third-Party (Optional)
- **`planet3D`** — renders 3D body meshes in orbit plots
  - https://www.mathworks.com/matlabcentral/fileexchange/86483
  - If not installed, comment out all `planet3D(...)` calls; orbits will still plot correctly
- **`violinplot`** — violin plot visualization for phase-point error distributions
  - https://www.mathworks.com/matlabcentral/fileexchange/45134

---

## Getting Started

### 1. Install dependencies
Install CVX and Gurobi and ensure both are on your MATLAB path. Verify with:
```matlab
cvx_setup
gurobi_setup
```

### 2. Configure the reference periodic orbit
At the top of `Fourier_Collocation_Template.m`, set the initial condition, period, and mass parameter for your target periodic orbit:

```matlab
ICp = [1.153555180492936, 0, 0.139438085327783, 0, -0.215370605225172, 0];
Tp  = 3.221357403048232;
u   = 0.012153593988908;
```

The provided values correspond to a vertical orbit near the Earth–Moon L2 point with Jacobi constant C ≈ 3.0776 LU²/TU².

### 3. Set algorithm parameters
Key tunable parameters are grouped near the top of each section:

```matlab
M            = 10;       % harmonic order for PO spectral decomposition
max_order_wP = 2;        % maximum order of quasi-periodic frequency harmonics
ds           = 5e-4;     % initial center manifold step size [LU]
t_LS         = 0:Tp/2:15*Tp;  % collocation arc endpoints
n_points     = 15;       % collocation sample points per arc
e_tol        = 0.0001;   % relative phase-point error tolerance for harmonic addition
index        = 4:4:36;   % QPO family indices to process in harmonic addition
```

### 4. Run
```matlab
run('Fourier_Collocation_Template.m')
```

The script is sequential. All outputs (figures, workspace variables) are generated in order. Expect the continuation loop and harmonic addition sections to be computationally intensive — runtimes of several hours are normal for larger families.

---

## Units

All states and times are in nondimensional CR3BP units unless explicitly converted for plotting:

| Quantity | Symbol | Value |
|---|---|---|
| Length unit | LU | 384,400 km (Earth–Moon distance) |
| Time unit | TU | 375,197 s ≈ 4.34 days |
| Frequency | rad/TU | — |

---

## Algorithm Summary

The full pipeline follows the flowchart in Fig. 6 of the paper:

```
Reference PO
    │
    ├─► Differential Correction  (enforce y(Tp)=0, ydot(Tp)=0)
    │
    ├─► Floquet Analysis          (extract ωₚ from monodromy matrix)
    │
    ├─► NAFF Decomposition        (fit Fourier series to PO trajectory)
    │
    ├─► Initial Torus             (seed QPO coefficients from PO fit)
    │
    ├─► Torus Refinement          (lsqnonlin: minimize ΔC, enforce collocation)
    │
    ├─► Arclength Continuation    (step through QPO family in Fourier space)
    │
    └─► Harmonic Addition         (NAFF on ΔC → CVX integer fit → re-solve)
              │
              └─► Phase-Point Drift Validation  (1-day propagation error)
```

### Continuation constraints
At each continuation step, three equality constraints enforce:
- **Arclength**: `(x - x_prev) · ξ = ds` — fixed step in Fourier space
- **Phase condition 1**: `(F - F_prev) · ∂q/∂Q₁ = 0` — no drift in Q₁ phase
- **Phase condition 2**: `(F - F_prev) · ∂q/∂Q₂ = 0` — no drift in Q₂ phase

And one inequality constraint enforces:
- **Fourier energy monotonicity**: `‖F_new‖² ≥ ‖F_prev‖²` — QPO amplitude must grow

### Error metrics (Eqs. 23–24 of the paper)
Phase-point drift is measured by launching a 50×50 grid of initial phases, propagating each for one day under true CR3BP dynamics, and comparing against the Fourier model:

```
E∞       = max_i ‖r̃_i − r_i‖₂                          (absolute position error)
(E∞)_rel = max_i ‖r̃_i − r_i‖₂ / ‖r_i − r_moon‖₂       (relative, normalized by Moon distance)
```

---

## Known Limitations

- **Eigenvector indexing** — The center manifold eigenvector is extracted at hardcoded indices `vec(:,2:3)` of the monodromy matrix. This is valid for the vertical orbit provided, but should be verified for other orbit families by inspecting the imaginary parts of `w_P` before proceeding.
- **Resonance handling** — The resonance detection threshold (`1e-3 rad/TU`) may require tuning for orbit families with different frequency scales.
- **Computational cost** — The collocation constraints launch `ode113` integrations inside `lsqnonlin`. For large harmonic bases or many arc segments, runtimes can be significant. Parallel evaluation (`UseParallel: true`) is enabled by default and strongly recommended.
- **Near-resonance behavior** — Continuation slows near rational frequency ratios ω₀/ωₚ ≈ integer. The algorithm reduces step size and retries, but may terminate early near strong resonances.

---

## Citation

Coming Soon

---

## Funding

This work was supported by the Department of Aerospace Engineering, Pennsylvania State University, and partially by the U.S. Air Force Office of Scientific Research (AFOSR) under grant FA9550-25-1-0078.

---

## License

This repository is released for academic and research use. Please contact the authors for other use cases.
