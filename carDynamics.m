function [dy, lambda, Sys_Input] = carDynamics(t, y, tau_func, mb, mw, Ibx, Iby, Ibz, Iw, g, R, l, w)
% Define symbolic variable ordering for q and qd (must match symbolic generation)
% vars_q  = {'x', 'y', 'z', 'phi', 'theta', 'psi', 'th1', 'th2', 'th3', 'th4'};
% vars_qd = {'xd', 'yd', 'zd', 'phid', 'thetad', 'psid', 'th1d', 'th2d', 'th3d', 'th4d'};

% Pack parameters into a vector matching matlabFunction order
params  = [mb; mw; Ibx; Iby; Ibz; Iw; g; R; l; w];

% Extract generalized coordinates and velocities from state vector
q  = y(1:12);
qd = y(13:24);


% Get input torques at current time
Tau = tau_func(t);

% Evaluate mass matrix M(q)
M = M_matrix_func(q, params);

% Evaluate Coriolis and centrifugal vector C(q, qd)
C = C_vector_func(q, qd, params);

% Evaluate gravity vector G(q)
G = G_vector_func(q, params);

% tau_gen = computeGeneralizedForces(q, Tau(7:10), params);  % computes forces + wheel torques
% qdd = M \ (tau_gen - C - G);

% ======= Augmented appraoch  ============
% Get Jacobians
[J,Jdot] = AllLegs_contactRolling_J_and_Jdot(q, qd, params);

[m,~] = size(J);


% Left-hand side (LHS)
LHS = [M, -J'; J, zeros(m,m)];
LHS = LHS + 1e-9 * eye(size(LHS));


% Right-hand side (RHS)
rhs = [Tau - C - G; -Jdot * qd];


% Solve for [qdd; lambda]
solution = LHS \ rhs;
qdd = solution(1:12);
lambda = solution(13:end);  % Contact forces lambda

Sys_Input = J'*lambda;

% Compose derivative of state vector
dy = [qd; qdd];
end


