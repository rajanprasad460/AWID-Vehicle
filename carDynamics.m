function [dy, lambda, Sys_Input] = carDynamics(t, y, tau_func, params)
% Define symbolic variable ordering for q and qd (must match symbolic generation)
% vars_q  = {'x', 'y', 'z', 'phi', 'theta', 'psi', 'th1', 'th2', 'th3', 'th4'};
% vars_qd = {'xd', 'yd', 'zd', 'phid', 'thetad', 'psid', 'th1d', 'th2d', 'th3d', 'th4d'};
% mb        = params(1);   % Mass of the body
% mw        = params(2);   % Mass of one wheel
% Ibx       = params(3);   % Body inertia about x-axis
% Iby       = params(4);   % Body inertia about y-axis
% Ibz       = params(5);   % Body inertia about z-axis
% Iw        = params(6);   % Inertia of a wheel
% g         = params(7);   % Gravitational acceleration
% R         = params(8);   % Radius of a wheel
% l         = params(9);   % Length (possibly between axles or wheelbase)
% w         = params(10);  % Width (track width or lateral offset)
% lat_damp  = params(11);  % Lateral damping coefficient




lat_damp  = params(11);  % Lateral damping coefficient
n_wheel = params(12);  % Number of wheels in the car


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
[J,Jdot,J_lat] = AllLegs_contactRolling_J_and_Jdot(q, qd, params);

% Initialize damping torque
tau_damping = (zeros(length(q),1));





for i = 1:n_wheel
    v_lat_i = J_lat(i,:) * qd;            % scalar lateral velocity of wheel i
    F_damp_i = -lat_damp * v_lat_i;         % damping force opposing lateral slip
    tau_damping = tau_damping + J_lat(i,:).' * F_damp_i; % add torque contribution
end



% tau_damping

[m,~] = size(J);



% Left-hand side (LHS)
LHS = [M, -J'; J, zeros(m,m)];
LHS = LHS + 1e-9 * eye(size(LHS));


% Right-hand side (RHS)
rhs = [(Tau - C - G -tau_damping) ; (-Jdot * qd)];


% Solve for [qdd; lambda]
solution = LHS \ rhs;
qdd = solution(1:12);
lambda = solution(13:end);  % Contact forces lambda

Sys_Input = J'*lambda;

% Compose derivative of state vector
dy = [qd; qdd];
end


