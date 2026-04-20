# Tubular Frame Structural Optimization

Optimization project for the minimum-mass design of a 2D tubular truss frame under stress, buckling, and geometric constraints. Implemented in MATLAB for the Optimization course (1st year, Q3).

## Problem Description

The project optimizes a **4-node, 5-element tubular truss** made of magnesium alloy. The goal is to minimize total structural mass while ensuring the frame remains safe under applied loads.

**Material:** Magnesium — density 0.0027 g/mm³, Young's modulus 69 GPa

**Applied loads:** 400 N at nodes 2–3 (vertical), −800 N at node 4

### Design Variables

| Problem | Variables |
|---|---|
| Simplified (2 vars) | Outer diameter D, inner diameter d |
| Complete (8 vars) | D, d, and x/y coordinates of nodes 2–4 |

### Constraints

- **Stress:** |σ| ≤ σ\_max
- **Buckling:** Euler critical load — P\_cr = π²EI / (KL)²
- **Geometry:** D − d ≥ 0.1 mm (wall thickness)
- **Strain:** ε ≤ 0.01%

## Methods Implemented

- **fmincon (SQP)** — gradient-based local optimization
- **Genetic Algorithm (GA)** — global search, 100 generations
- **Hybrid GA + fmincon** — global search refined locally
- **Multistart** — 6×6 grid of initial points to map convergence basins
- **Steepest Descent** — manual implementation with finite-difference gradients and line search
- **Finite Difference Analysis** — gradient accuracy vs. step size study

## Repository Structure

```
simplified_problem.m    # 2-variable optimization (D, d only)
complete_problem.m      # 8-variable optimization (D, d + node positions)
final_optimization-6298745.pdf  # Full report with analysis and results
```

## Structural Analysis

The frame is solved using 2D truss equilibrium (A·T = F). Cross-sectional properties:

- Area: S = π(D² − d²) / 4
- Moment of inertia: I = (π/64)(D⁴ − d⁴)

The penalty method (exponential) is used to handle constraint violations during unconstrained search phases.

## Outputs

Each optimization method reports:
- Optimal mass and design variables
- Internal forces and stresses in each element
- Feasibility of all constraints
- Comparison across methods

Visualizations include 3D surface and contour plots of the objective function, constraint boundaries, multistart convergence scatter plots, and deformed frame geometry.

## Requirements

- MATLAB R2020b or later
- Optimization Toolbox (for `fmincon`, `ga`)
- Global Optimization Toolbox (for `MultiStart`)
