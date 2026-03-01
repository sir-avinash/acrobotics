# Acrobotics

Dynamics simulation and comparison of a **two-link acrobot** (acrobatic robot / double pendulum) using three different mathematical representations of rigid body orientation:

1. **Euler Angles** — Classical roll/pitch parametrization
2. **S² (Unit Sphere)** — Direction vectors on the 2-sphere with angular velocity
3. **SO(3) (Special Orthogonal Group)** — Full rotation matrices with body-frame angular velocity

The project derives equations of motion using Lagrangian mechanics and integrates them with MATLAB's `ode45` solver, then compares kinetic energy, potential energy, angles, direction vectors, and angular velocities across all three representations.

## System Model

The acrobot consists of two rigid links connected in series:

- **Link 1**: mass `m1`, length `l1`, hanging from a fixed pivot
- **Link 2**: mass `m2`, length `l2`, connected to the end of link 1
- Both links are modeled as **point masses** at the link endpoints (Euler and S² models)
- The SO(3) model additionally includes a small rotational inertia `J = 1e-3 * diag([0,0,1])`
- Gravity acts along the **+e3** direction (z-axis up); the pendulum starts nearly vertical

The system is **unactuated** (no control input) — it evolves under gravity alone.

### Default Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `g`       | 9.81  | Gravitational acceleration (m/s²) |
| `m`       | 1     | Mass per link (kg) |
| `l`       | 1     | Link length (m) |
| `J`       | 1e-3·diag([0,0,1]) | Moment of inertia (SO(3) only) |
| `theta`   | -10°  | Initial pitch (both links) |
| `phi`     | 0°    | Initial roll (both links) |
| `Tend`    | 1 s   | Simulation duration |
| `Fs`      | 80 Hz | Sampling frequency for plots |

## File Structure

| File | Description |
|------|-------------|
| `acrobotics_main.m` | Main simulation script — sets up parameters, runs ODE integration for all three representations, computes energies, and generates comparison plots |
| `derive_acrobot_eulerdynamics.m` | Symbolic derivation of Euler-angle equations of motion using Lagrangian mechanics (requires Symbolic Math Toolbox) |
| `acrobot_inertia_matrix.m` | Auto-generated 4×4 mass/inertia matrix `D(q)` for Euler representation |
| `acrobot_coriolis_and_gravity.m` | Auto-generated Coriolis + gravity vector `H(q, dq)` for Euler representation |
| `acrobot_kinetic_energy.m` | Auto-generated kinetic energy function `KE(q, dq)` |
| `acrobot_potential_energy.m` | Auto-generated potential energy function `PE(q)` |
| `vecnorm.m` | Utility: column-wise vector norms of a matrix |

## Dependencies

- **MATLAB** (tested with Symbolic Math Toolbox v7.2, October 2017)
- **External Toolbox** providing: `hat` (skew-symmetric matrix), `Rx`/`Ry` (rotation matrices), `Rot2Eul` (rotation to Euler angles), `cross2` (cross product), `even_sample` (uniform time resampling), `EoM` (Euler-Lagrange equation derivation)

> **Note**: The `addpath` in `acrobotics_main.m` (line 5) points to a machine-specific path and must be updated to the location of the external toolbox on your system.

## Usage

### Running the Simulation

```matlab
% 1. Update the toolbox path in acrobotics_main.m (line 5)
% 2. Run:
acrobotics_main
```

This produces 6 comparison figures:
- **Figure 1**: Kinetic, potential, and total energy for all three representations
- **Figure 2**: Roll (φ) and pitch (θ) angles for both links
- **Figures 4–5**: Direction vector components (q1, q2) and their norms
- **Figures 6–7**: Angular velocity components (ω1, ω2)

### Regenerating Euler Dynamics

To re-derive the symbolic equations of motion and regenerate the `acrobot_*.m` functions:

```matlab
derive_acrobot_eulerdynamics
```

This overwrites `acrobot_inertia_matrix.m`, `acrobot_coriolis_and_gravity.m`, `acrobot_kinetic_energy.m`, and `acrobot_potential_energy.m`.

## Mathematical Formulation

### Euler Angles

Generalized coordinates: `q = [φ₁, θ₁, φ₂, θ₂]` (roll and pitch for each link).
Direction vectors are recovered via `qᵢ = Rx(φᵢ)·Ry(θᵢ)·e₃`.
Equations of motion: `D(q)·q̈ = H(q, q̇) + B·u`, where `H = -C·q̇ - G` lumps Coriolis and gravity terms.

### S² (Unit Sphere)

State: `(q₁, q₂, Ω₁, Ω₂)` where `qᵢ ∈ S²` are unit direction vectors and `Ωᵢ` are angular velocities satisfying `q̇ᵢ = Ωᵢ × qᵢ`.
The dynamics are formulated as constrained Euler-Lagrange equations on the sphere.

### SO(3) (Rotation Group)

State: `(R₁, R₂, ω₁, ω₂)` where `Rᵢ ∈ SO(3)` are rotation matrices and `ωᵢ` are angular velocities.
This representation includes rotational inertia `J` and avoids singularities present in Euler angles.

## License

Licensed under the [Apache License, Version 2.0](LICENSE).
