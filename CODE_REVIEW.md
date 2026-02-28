# Code Review: Acrobotics

## Critical Issues

### 1. SO(3) Rotation Kinematics Bug (`acrobotics_main.m:481-482`)

The rotation matrix derivative uses the **wrong formula** for a body-to-world rotation matrix.

```matlab
% Current (WRONG):
R1_dot = -hat(omega1)*R1w;   % line 481
R2_dot = -hat(omega2)*R2w;   % line 482
```

From the initial conditions, `R1w = Rx(phi)*Ry(theta)` maps body frame to world frame (since `q = R1w*e3` gives the pendulum direction in world coordinates). For a body-to-world rotation matrix, the correct kinematic equation is:

```matlab
% Correct:
R1_dot = R1w * hat(omega1);   % dR_b2w/dt = R_b2w * hat(omega_body)
R2_dot = R2w * hat(omega2);
```

The formula `dR/dt = -hat(omega_body)*R` is only valid for a **world-to-body** rotation matrix. Using it on a body-to-world matrix produces `dq/dt = -(omega_body x q)` instead of the correct `dq/dt = omega_world x q`.

### 2. SO(3) Angular Acceleration Frame Transformation Bug (`acrobotics_main.m:484-485`)

The angular acceleration is transformed with the **inverse** of the correct rotation.

```matlab
% Current (WRONG):
omega1_dot = R1*out(1:3);    % R1 = R1w' (world-to-body)
omega2_dot = R2*out(4:6);
```

Since the state stores `omega_world` and the dynamics solve for `alpha_body`:
- `omega_world_dot = R_b2w * alpha_body` (since `omega_body x omega_body = 0`)
- This requires `R1w`, not `R1 = R1w'`

```matlab
% Correct:
omega1_dot = R1w * out(1:3);   % R_b2w * alpha_body
omega2_dot = R2w * out(4:6);
```

**Impact**: These two bugs cause the SO(3) simulation to produce incorrect trajectories. For small angles and short time horizons (the default `Tend=1`), the errors may be small enough to appear reasonable on plots, but they will grow over time.

### 3. Data Save Logic Never Triggers (`acrobotics_main.m:195-197`)

```matlab
if exist('main_data')==1 && force_save==1
    save('main_data')
end
```

`exist('main_data')==1` checks whether `main_data` is a **variable** in the workspace. No such variable is ever created, so this condition is always false and the data is never saved, even though `force_save=1`. The logical operator should likely be `||` (save if forced OR if data already exists), or the condition should simply be:

```matlab
if force_save == 1
    save('main_data')
end
```

---

## Moderate Issues

### 4. Euler Angle Recovery Singularities (`acrobotics_main.m:104-107`)

The formulas for recovering `dth` and `dph` from S^2 states involve divisions that blow up at certain angles:

```matlab
spherical.dth1(i) = spherical.dq1(1,i)/cos(spherical.th1(i));          % singular at theta = +/-pi/2
spherical.dph1(i) = ... / (cos(th)*(sin(ph)+cos(ph)));                  % singular at phi = -pi/4 + n*pi
```

This will produce `Inf`/`NaN` values if the pendulum swings through these configurations. A more robust approach would use the full `dq` vector with a pseudoinverse or quaternion-based recovery.

### 5. Unused Variable `mt` in Euler Simulation (`acrobotics_main.m:406`)

```matlab
mt = data.mt;   % extracted but never used in acrobot_sim_euler
```

The torso mass `mt` is extracted from the data struct but plays no role in the Euler dynamics function. It appears to be a leftover from a previous version (the `.asv` backup has `data.mt = 2`). This dead code is misleading about what the Euler model actually simulates.

### 6. SO(3) Dynamics Include Rotational Inertia, Others Don't

The SO(3) model uses `data.J = 1e-3*diag([0,0,1])` for rotational inertia of each link, while the S^2 and Euler models treat links as **point masses** (no rotational inertia). This means the three representations are modeling **slightly different physical systems**, which undermines the comparison. With `J = 1e-3*diag([0,0,1])`, the effect is small but non-zero.

### 7. `vecnorm.m` Shadows MATLAB Built-in (`vecnorm.m`)

Starting from MATLAB R2017b, `vecnorm` is a built-in function. This custom implementation shadows the built-in, which can cause subtle issues if other toolbox code expects the built-in behavior (different calling convention and output shape).

---

## Minor Issues / Code Quality

### 8. Hardcoded Toolbox Path (`acrobotics_main.m:5`)

```matlab
addpath('/home/exx/Avinash/My Toolbox');
```

This absolute path is machine-specific and will fail on any other system. The external toolbox dependency (`hat`, `Rx`, `Ry`, `Rot2Eul`, `cross2`, `even_sample`, `EoM`) should be documented and the path made configurable.

### 9. Figure 3 is Skipped

Plots use figures 1, 2, 4, 5, 6, 7 but skip figure 3. This appears to be from commented-out code that previously used figure 3.

### 10. Backup File in Repository (`acrobotics_main.asv`)

MATLAB auto-save `.asv` files should not be committed. This file contains a stale version with different parameters (`data.mt = 2`, `data.J = 0*eye(3)`, `Tend = 4`).

### 11. Empty/Stub Files

- `acrobot_input.m` - empty (0 bytes)
- `main_s1` - empty
- `main_so2` - empty
- `main_euler.m` - only a comment

These appear to be placeholders for future work that was never completed.

### 12. Large Amounts of Commented-Out Code

Lines 282-367 of `acrobotics_main.m` contain ~85 lines of commented-out plotting code. This should be removed or moved to a separate analysis script.

### 13. `clear all` is Commented Out (`acrobotics_main.m:3`)

```matlab
% clear all
```

Running the script multiple times without clearing the workspace can lead to stale variable contamination, especially since the save condition (issue #3) checks for workspace variables.

---

## Summary

| Severity | Count | Issues |
|----------|-------|--------|
| Critical | 2 | SO(3) kinematics bugs (#1, #2) |
| Moderate | 5 | Save logic (#3), singularities (#4), unused var (#5), model mismatch (#6), name shadowing (#7) |
| Minor    | 6 | Hardcoded path (#8), skipped figure (#9), .asv in repo (#10), stubs (#11), dead code (#12-13) |

The most impactful bugs are in the SO(3) simulation function where the rotation matrix derivative and angular acceleration transformation use incorrect frame conventions. These cause the SO(3) representation to diverge from the physically correct solution, especially for larger angles or longer simulation times.
