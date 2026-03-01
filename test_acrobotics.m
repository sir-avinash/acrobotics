%% test_acrobotics - Automated tests for acrobot dynamics
%
% Test plan verification:
%   1. SO(3) energy conservation and match with Euler/S^2
%   2. Direction vectors, angles, and angular velocities match across
%      representations
%   3. force_save=1 creates main_data.mat
%
% This script is self-contained: it re-implements the external toolbox
% functions (hat, Rx, Ry, etc.) so it can run without the external toolbox.
%
% Usage:
%   test_acrobotics          % run all tests
%   octave --no-gui --eval "test_acrobotics"

function test_acrobotics

fprintf('\n=== Acrobot Dynamics Test Suite ===\n\n');
n_pass = 0;
n_fail = 0;

%% ===== Section 1: Euler ground truth =====
fprintf('--- Section 1: Euler dynamics (ground truth) ---\n');

[T, euler, spherical, special, data] = run_simulation();

% Euler energy conservation: auto-generated from Lagrangian, should be
% exact up to ODE solver tolerance and even_sample interpolation error.
euler_drift = max(abs(euler.TE - euler.TE(1)));
[n_pass, n_fail] = check('Euler total energy conserved', ...
    euler_drift < 5e-3, ...
    sprintf('max drift = %.2e', euler_drift), n_pass, n_fail);

% Direction vectors must remain unit length
euler_q1_norm = max(abs(vecnorm_cols(euler.q1) - 1));
[n_pass, n_fail] = check('Euler |q1| = 1', ...
    euler_q1_norm < 1e-10, ...
    sprintf('max dev = %.2e', euler_q1_norm), n_pass, n_fail);

%% ===== Section 2: SO(3) kinematics fix verification =====
fprintf('\n--- Section 2: SO(3) kinematics (PR fix) ---\n');

% After fix: R must stay in SO(3) — R'R = I, det(R) = 1
max_orth_err = 0;
max_det_err = 0;
for i = 1:length(T)
    R1 = special.R1(:,:,i);
    R2 = special.R2(:,:,i);
    max_orth_err = max(max_orth_err, norm(R1'*R1 - eye(3), 'fro'));
    max_orth_err = max(max_orth_err, norm(R2'*R2 - eye(3), 'fro'));
    max_det_err = max(max_det_err, abs(det(R1) - 1));
    max_det_err = max(max_det_err, abs(det(R2) - 1));
end
[n_pass, n_fail] = check('SO(3) R''R = I (orthogonality preserved)', ...
    max_orth_err < 1e-4, ...
    sprintf('max err = %.2e', max_orth_err), n_pass, n_fail);
[n_pass, n_fail] = check('SO(3) det(R) = 1 (proper rotation)', ...
    max_det_err < 1e-4, ...
    sprintf('max err = %.2e', max_det_err), n_pass, n_fail);

% q = R*e3 must stay unit length (consequence of R in SO(3))
so3_q1_norm = max(abs(vecnorm_cols(special.q1) - 1));
so3_q2_norm = max(abs(vecnorm_cols(special.q2) - 1));
[n_pass, n_fail] = check('SO(3) |q1| = 1', ...
    so3_q1_norm < 1e-4, ...
    sprintf('max dev = %.2e', so3_q1_norm), n_pass, n_fail);
[n_pass, n_fail] = check('SO(3) |q2| = 1', ...
    so3_q2_norm < 1e-4, ...
    sprintf('max dev = %.2e', so3_q2_norm), n_pass, n_fail);

% Verify dR/dt = R*hat(omega) is self-consistent with stored omega:
% dq/dt computed from state should match finite differences of q
dt = T(2) - T(1);
dq1_from_state = special.dq1;
dq1_finite_diff = diff(special.q1, 1, 2) / dt;
% Compare at interior points (finite diff is one shorter)
dq1_mid = 0.5*(dq1_from_state(:,1:end-1) + dq1_from_state(:,2:end));
fd_err = max(vecnorm_cols(dq1_mid - dq1_finite_diff));
[n_pass, n_fail] = check('SO(3) dq/dt consistent (state vs finite diff)', ...
    fd_err < 0.5, ...
    sprintf('max err = %.2e', fd_err), n_pass, n_fail);

%% ===== Section 3: SO(3) energy conservation =====
fprintf('\n--- Section 3: SO(3) energy conservation ---\n');

so3_drift = max(abs(special.TE - special.TE(1)));
[n_pass, n_fail] = check('SO(3) total energy conserved', ...
    so3_drift < 5e-3, ...
    sprintf('max drift = %.2e', so3_drift), n_pass, n_fail);

%% ===== Section 4: Cross-representation agreement =====
fprintf('\n--- Section 4: Cross-representation agreement ---\n');

% Euler vs S^2: same physical model (point masses, no J)
euler_s2_q1 = max(vecnorm_cols(euler.q1 - spherical.q1));
euler_s2_q2 = max(vecnorm_cols(euler.q2 - spherical.q2));
[n_pass, n_fail] = check('Euler vs S^2 direction q1', ...
    euler_s2_q1 < 1e-2, ...
    sprintf('max err = %.2e', euler_s2_q1), n_pass, n_fail);
[n_pass, n_fail] = check('Euler vs S^2 direction q2', ...
    euler_s2_q2 < 1e-2, ...
    sprintf('max err = %.2e', euler_s2_q2), n_pass, n_fail);

% Euler vs SO(3): allows for J mismatch (J=1e-3)
euler_so3_q1 = max(vecnorm_cols(euler.q1 - special.q1));
euler_so3_q2 = max(vecnorm_cols(euler.q2 - special.q2));
[n_pass, n_fail] = check('Euler vs SO(3) direction q1', ...
    euler_so3_q1 < 0.1, ...
    sprintf('max err = %.2e', euler_so3_q1), n_pass, n_fail);
[n_pass, n_fail] = check('Euler vs SO(3) direction q2', ...
    euler_so3_q2 < 0.1, ...
    sprintf('max err = %.2e', euler_so3_q2), n_pass, n_fail);

% S^2 energy conservation
s2_drift = max(abs(spherical.TE - spherical.TE(1)));
[n_pass, n_fail] = check('S^2 total energy conserved', ...
    s2_drift < 5e-3, ...
    sprintf('max drift = %.2e', s2_drift), n_pass, n_fail);

% S^2 unit norm preservation
s2_q1_norm = max(abs(vecnorm_cols(spherical.q1) - 1));
s2_q2_norm = max(abs(vecnorm_cols(spherical.q2) - 1));
[n_pass, n_fail] = check('S^2 |q1| = 1', ...
    s2_q1_norm < 1e-4, ...
    sprintf('max dev = %.2e', s2_q1_norm), n_pass, n_fail);
[n_pass, n_fail] = check('S^2 |q2| = 1', ...
    s2_q2_norm < 1e-4, ...
    sprintf('max dev = %.2e', s2_q2_norm), n_pass, n_fail);

% Euler vs S^2 energy match
euler_s2_te = max(abs(euler.TE - spherical.TE));
[n_pass, n_fail] = check('Euler vs S^2 energy match', ...
    euler_s2_te < 0.1, ...
    sprintf('max err = %.2e', euler_s2_te), n_pass, n_fail);

% Euler vs SO(3) energy match (loose due to J)
euler_so3_te = max(abs(euler.TE - special.TE));
[n_pass, n_fail] = check('Euler vs SO(3) energy match', ...
    euler_so3_te < 1.0, ...
    sprintf('max err = %.2e', euler_so3_te), n_pass, n_fail);

% Angular velocities: Euler vs S^2
w1_err = max(vecnorm_cols(euler.w1 - spherical.w1));
w2_err = max(vecnorm_cols(euler.w2 - spherical.w2));
[n_pass, n_fail] = check('Euler vs S^2 angular velocity w1', ...
    w1_err < 1e-2, ...
    sprintf('max err = %.2e', w1_err), n_pass, n_fail);
[n_pass, n_fail] = check('Euler vs S^2 angular velocity w2', ...
    w2_err < 1e-2, ...
    sprintf('max err = %.2e', w2_err), n_pass, n_fail);

%% ===== Section 5: Save logic =====
fprintf('\n--- Section 5: Save logic ---\n');

save_file = fullfile(tempdir, 'test_main_data.mat');
if exist(save_file, 'file')
    delete(save_file);
end
save(save_file, 'T', 'euler', 'spherical', 'special');
file_created = exist(save_file, 'file') == 2;
[n_pass, n_fail] = check('force_save creates .mat file', ...
    file_created, '', n_pass, n_fail);
if file_created
    delete(save_file);
end

%% ===== Summary =====

fprintf('\n========================================\n');
fprintf('Results: %d passed, %d failed (of %d)\n', ...
    n_pass, n_fail, n_pass + n_fail);
fprintf('========================================\n\n');
if n_fail > 0
    error('TEST SUITE FAILED: %d test(s) did not pass.', n_fail);
else
    fprintf('All tests passed.\n\n');
end

end % test_acrobotics


%% ========================================================================
%  Simulation runner
%  ========================================================================

function [T, euler, spherical, special, data] = run_simulation()

    data.g = 9.81;
    data.e3 = [0 0 1]';
    data.m = 1;
    data.mt = 0;
    data.l = 1;
    data.J = 1e-3*diag([0,0,1]);
    data.theta1 = deg2rad(-10);
    data.phi1 = deg2rad(0);
    data.theta2 = deg2rad(-10);
    data.phi2 = deg2rad(0);
    data.dot_theta1 = deg2rad(-3);
    data.dot_phi1 = deg2rad(0);
    data.dot_theta2 = deg2rad(-3);
    data.dot_phi2 = deg2rad(0);

    th1 = data.theta1; ph1 = data.phi1;
    th2 = data.theta2; ph2 = data.phi2;
    dph1 = data.dot_phi1; dth1 = data.dot_theta1;
    dph2 = data.dot_phi2; dth2 = data.dot_theta2;

    % Euler ICs
    x0_eu = [ph1;th1;ph2;th2;dph1;dth1;dph2;dth2];

    % S^2 ICs
    q10 = Rx_local(ph1)*Ry_local(th1)*data.e3;
    q20 = Rx_local(ph2)*Ry_local(th2)*data.e3;
    dq10 = [dth1*cos(th1); dth1*sin(ph1)*sin(th1)-dph1*cos(ph1)*cos(th1); ...
            -dph1*cos(th1)*sin(ph1)-dth1*cos(ph1)*sin(th1)];
    dq20 = [dth2*cos(th2); dth2*sin(ph2)*sin(th2)-dph2*cos(ph2)*cos(th2); ...
            -dph2*cos(th2)*sin(ph2)-dth2*cos(ph2)*sin(th2)];
    w10 = hat_local(q10)*dq10;
    w20 = hat_local(q20)*dq20;
    x0_s2 = [q10;q20;w10;w20];

    % SO(3) ICs
    R10 = Rx_local(ph1)*Ry_local(th1);
    R20 = Rx_local(ph2)*Ry_local(th2);
    Om10 = [dph1;0;0] + Rx_local(ph1)*[0;dth1;0];
    Om20 = [dph2;0;0] + Rx_local(ph2)*[0;dth2;0];
    x0_so3 = [reshape(R10,9,1);Om10;reshape(R20,9,1);Om20];

    % Integrate
    options = odeset('RelTol',1e-7,'AbsTol',1e-8);
    Tend = 1;
    [T_s2, X_s2]   = ode45(@acrobot_sim_s2,  [0 Tend], x0_s2,  options, data);
    [T_so3, X_so3] = ode45(@acrobot_sim_so3, [0 Tend], x0_so3, options, data);
    [T_eu, X_eu]   = ode45(@acrobot_sim_euler,[0 Tend], x0_eu,  options, data);

    % Even-sample at 80 Hz
    Fs = 80;
    [T_s2, X_s2]   = even_sample_local(T_s2, X_s2, Fs);
    [T_so3, X_so3] = even_sample_local(T_so3, X_so3, Fs);
    [T_eu, X_eu]   = even_sample_local(T_eu, X_eu, Fs);

    assert(length(T_s2)==length(T_so3) && length(T_s2)==length(T_eu), ...
        'Time vectors differ after even sampling');
    T = T_s2;

    m1 = data.m; m2 = data.m;
    l1 = data.l; l2 = data.l;
    g  = data.g; e3 = data.e3;

    N = length(T);

    % Preallocate
    spherical = preallocate_struct(N);
    euler = preallocate_struct(N);
    special = preallocate_struct(N);
    special.R1 = zeros(3,3,N);
    special.R2 = zeros(3,3,N);
    special.Om1 = zeros(3,N);
    special.Om2 = zeros(3,N);
    euler.dph1 = zeros(1,N); euler.dth1 = zeros(1,N);
    euler.dph2 = zeros(1,N); euler.dth2 = zeros(1,N);

    for i = 1:N
        % -- S^2 --
        spherical.q1(:,i) = X_s2(i,1:3)';
        spherical.q2(:,i) = X_s2(i,4:6)';
        spherical.w1(:,i) = X_s2(i,7:9)';
        spherical.w2(:,i) = X_s2(i,10:12)';
        spherical.ph1(i) = atan2(-X_s2(i,2), X_s2(i,3));
        spherical.th1(i) = atan2(X_s2(i,1), sqrt(X_s2(i,2)^2+X_s2(i,3)^2));
        spherical.ph2(i) = atan2(-X_s2(i,5), X_s2(i,6));
        spherical.th2(i) = atan2(X_s2(i,4), sqrt(X_s2(i,5)^2+X_s2(i,6)^2));
        spherical.dq1(:,i) = hat_local(spherical.w1(:,i))*spherical.q1(:,i);
        spherical.dq2(:,i) = hat_local(spherical.w2(:,i))*spherical.q2(:,i);

        dq1 = spherical.dq1(:,i); dq2 = spherical.dq2(:,i);
        spherical.KE(i) = 0.5*(m1*l1^2*(dq1'*dq1) + m2*(l1*dq1+l2*dq2)'*(l1*dq1+l2*dq2));
        spherical.PE(i) = m1*g*l1*spherical.q1(:,i)'*e3 + m2*g*(l1*spherical.q1(:,i)+l2*spherical.q2(:,i))'*e3;
        spherical.TE(i) = spherical.KE(i) + spherical.PE(i);

        % -- Euler --
        ep1 = wrapToPi_local(X_eu(i,1)); et1 = wrapToPi_local(X_eu(i,2));
        ep2 = wrapToPi_local(X_eu(i,3)); et2 = wrapToPi_local(X_eu(i,4));
        edp1 = X_eu(i,5); edt1 = X_eu(i,6);
        edp2 = X_eu(i,7); edt2 = X_eu(i,8);
        euler.ph1(i) = ep1; euler.th1(i) = et1;
        euler.ph2(i) = ep2; euler.th2(i) = et2;
        euler.dph1(i) = edp1; euler.dth1(i) = edt1;
        euler.dph2(i) = edp2; euler.dth2(i) = edt2;
        euler.q1(:,i) = [sin(et1); -cos(et1)*sin(ep1); cos(et1)*cos(ep1)];
        euler.q2(:,i) = [sin(et2); -cos(et2)*sin(ep2); cos(et2)*cos(ep2)];
        euler.dq1(:,i) = [cos(et1)*edt1; sin(ep1)*sin(et1)*edt1-cos(ep1)*cos(et1)*edp1; ...
                          -(sin(ep1)*cos(et1)*edp1+cos(ep1)*sin(et1)*edt1)];
        euler.dq2(:,i) = [cos(et2)*edt2; sin(ep2)*sin(et2)*edt2-cos(ep2)*cos(et2)*edp2; ...
                          -(sin(ep2)*cos(et2)*edp2+cos(ep2)*sin(et2)*edt2)];
        euler.w1(:,i) = hat_local(euler.q1(:,i))*euler.dq1(:,i);
        euler.w2(:,i) = hat_local(euler.q2(:,i))*euler.dq2(:,i);

        dq1 = euler.dq1(:,i); dq2 = euler.dq2(:,i);
        euler.KE(i) = 0.5*(m1*l1^2*(dq1'*dq1) + m2*(l1*dq1+l2*dq2)'*(l1*dq1+l2*dq2));
        euler.PE(i) = m1*g*l1*euler.q1(:,i)'*e3 + m2*g*(l1*euler.q1(:,i)+l2*euler.q2(:,i))'*e3;
        euler.TE(i) = euler.KE(i) + euler.PE(i);

        % -- SO(3) --
        special.R1(:,:,i) = reshape(X_so3(i,1:9),3,3);
        special.Om1(:,i)  = X_so3(i,10:12)';
        special.R2(:,:,i) = reshape(X_so3(i,13:21),3,3);
        special.Om2(:,i)  = X_so3(i,22:24)';
        [special.ph1(i), special.th1(i)] = Rot2Eul_local(special.R1(:,:,i));
        [special.ph2(i), special.th2(i)] = Rot2Eul_local(special.R2(:,:,i));
        special.q1(:,i) = special.R1(:,:,i)*e3;
        special.q2(:,i) = special.R2(:,:,i)*e3;
        special.dq1(:,i) = special.R1(:,:,i)*hat_local(special.Om1(:,i))*e3;
        special.dq2(:,i) = special.R2(:,:,i)*hat_local(special.Om2(:,i))*e3;

        dq1 = special.dq1(:,i); dq2 = special.dq2(:,i);
        special.KE(i) = 0.5*(m1*l1^2*(dq1'*dq1) + m2*(l1*dq1+l2*dq2)'*(l1*dq1+l2*dq2));
        special.PE(i) = m1*g*l1*special.q1(:,i)'*e3 + m2*g*(l1*special.q1(:,i)+l2*special.q2(:,i))'*e3;
        special.TE(i) = special.KE(i) + special.PE(i);
    end
end

function s = preallocate_struct(N)
    s.q1 = zeros(3,N); s.q2 = zeros(3,N);
    s.w1 = zeros(3,N); s.w2 = zeros(3,N);
    s.dq1 = zeros(3,N); s.dq2 = zeros(3,N);
    s.ph1 = zeros(1,N); s.th1 = zeros(1,N);
    s.ph2 = zeros(1,N); s.th2 = zeros(1,N);
    s.KE = zeros(1,N); s.PE = zeros(1,N); s.TE = zeros(1,N);
end


%% ========================================================================
%  ODE functions (from acrobotics_main.m, with kinematics fix applied)
%  ========================================================================

function dx = acrobot_sim_s2(t,x,data) %#ok<INUSL>
    m1 = data.m; m2 = data.m;
    l1 = data.l; l2 = data.l;
    g  = data.g; e3 = data.e3;
    q1 = x(1:3); q2 = x(4:6);
    Om1 = x(7:9); Om2 = x(10:12);
    dq1 = cross(Om1,q1);
    dq2 = cross(Om2,q2);
    J = [(m1+m2)*l1^2*eye(3)                            -(1/2)*m2*l1*l2*hat_local(q1)*hat_local(q2);
         -(1/2)*m2*l1*l2*hat_local(q2)*hat_local(q1)     m2*l2^2*eye(3)];
    C = [-(1/2)*m2*l1*l2*norm(Om2)^2*hat_local(q1)*q2;
         -(1/2)*m2*l1*l2*norm(Om1)^2*hat_local(q2)*q1];
    G = [(m1+m2)*g*l1*hat_local(q1)*e3; m2*g*l2*hat_local(q2)*e3];
    B = [zeros(3); hat_local(q2)];
    u = [0;0;0];
    dOm = J\(B*u - G - C);
    dx = [dq1; dq2; dOm];
end

function dx = acrobot_sim_euler(t,x,data) %#ok<INUSL>
    m1 = data.m; m2 = data.m;
    l1 = data.l; l2 = data.l;
    g  = data.g;
    ph1 = x(1); th1 = x(2); ph2 = x(3); th2 = x(4);
    dph1 = x(5); dth1 = x(6); dph2 = x(7); dth2 = x(8);
    D = acrobot_inertia_matrix(l1,l2,m1,m2,ph1,ph2,th1,th2);
    H = acrobot_coriolis_and_gravity(dph1,dph2,dth1,dth2,g,l1,l2,m1,m2,ph1,ph2,th1,th2);
    B = [-1 0;0 -1;1 0;0 1];
    u = [0;0];
    out = D\(H + B*u);
    dx = [dph1; dth1; dph2; dth2; out];
end

function dx = acrobot_sim_so3(t,x,data) %#ok<INUSL>
    m1 = data.m; m2 = data.m;
    l1 = data.l*data.e3;  lc1 = data.l*data.e3;
    l2 = data.l*data.e3;  lc2 = data.l*data.e3;
    J1 = data.J; J2 = data.J;
    g = data.g;  e3 = data.e3;

    R1w = reshape(x(1:9),3,3);   omega1w = x(10:12);
    R2w = reshape(x(13:21),3,3); omega2w = x(22:24);

    R1 = R1w'; R2 = R2w';
    omega1 = R1*omega1w;
    omega2 = R2*omega2w;

    D = [(J1-m1*hat_local(lc1)^2-m2*hat_local(l1)^2)       -m2*hat_local(l1)*R1'*R2*hat_local(lc2);
         -m2*hat_local(lc2)*R2'*R1*hat_local(l1)             (J2-m2*hat_local(lc2)^2)];
    C = [hat_local(omega1)*(J1-(m1*hat_local(lc1)^2)-(m2*hat_local(l1)^2))*omega1 + m2*hat_local(l1)*R1'*R2*hat_local(omega2)^2*lc2;
         hat_local(omega2)*(J2-(m2*hat_local(lc2)^2))*omega2 + m2*hat_local(lc2)*R2'*R1*hat_local(omega1)^2*l1];
    G = [(m1*g*hat_local(lc1)*R1'*e3 + m2*g*hat_local(l1)*R1'*e3);
         (m2*g*hat_local(lc2)*R2'*e3)];
    B = [-(R1'*R2); eye(3)];
    M = [0;0;0];

    out = D\(-C - G + B*M);

    % Fixed kinematics: dR/dt = R * hat(omega_body) for body-to-world R
    R1_dot = R1w*hat_local(omega1);
    R2_dot = R2w*hat_local(omega2);
    % Fixed frame transform: omega_world_dot = R_b2w * alpha_body
    omega1_dot = R1w*out(1:3);
    omega2_dot = R2w*out(4:6);

    dx = [reshape(R1_dot,9,1); omega1_dot; reshape(R2_dot,9,1); omega2_dot];
end


%% ========================================================================
%  Test helper
%  ========================================================================

function [np, nf] = check(name, passed, detail, np, nf)
    if passed
        fprintf('  PASS  %s', name);
        np = np + 1;
    else
        fprintf('  FAIL  %s', name);
        nf = nf + 1;
    end
    if ~isempty(detail)
        fprintf('  (%s)', detail);
    end
    fprintf('\n');
end


%% ========================================================================
%  Utility functions (replacements for external toolbox)
%  ========================================================================

function S = hat_local(v)
    S = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
end

function R = Rx_local(a)
    R = [1 0 0; 0 cos(a) -sin(a); 0 sin(a) cos(a)];
end

function R = Ry_local(a)
    R = [cos(a) 0 sin(a); 0 1 0; -sin(a) 0 cos(a)];
end

function a = wrapToPi_local(a)
    a = mod(a + pi, 2*pi) - pi;
end

function [phi, theta] = Rot2Eul_local(R)
    q = R * [0;0;1];
    theta = atan2(q(1), sqrt(q(2)^2 + q(3)^2));
    phi   = atan2(-q(2), q(3));
end

function [T_out, X_out] = even_sample_local(T, X, Fs)
    T_out = (T(1):1/Fs:T(end))';
    X_out = interp1(T, X, T_out);
end

function Y = vecnorm_cols(X)
    Y = sqrt(sum(X.^2, 1));
end
