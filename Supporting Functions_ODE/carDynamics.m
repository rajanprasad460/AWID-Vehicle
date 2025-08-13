function dy = carDynamics(t, y, p, tau_func, params, reset)
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
% roll_damp  = params(12);  % rolling damping coefficient
% n_wheel = params(13);  % Number of wheels in the car




% lat_damp  = params(11);  % Lateral damping coefficient
% roll_damp  = params(12);  % rolling damping coefficient
% n_wheel = params(13);  % Number of wheels in the car



% === Persistent variables ===
% Used for timing logs, energy tracking, and cable force history
persistent lastTime

% === Reset persistent state ===
if nargin > 5 && reset
    lastTime = [];
    return
end
% === Initialize persistent variables ===
if isempty(lastTime)
    lastTime = tic;   % Start timer
end

% Extract generalized coordinates and velocities from state vector
q  = y(1:12);
qd = y(13:24);


% Get input torques at current time
Tau = tau_func(t, y);

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
[J,Jdot,J_Total,~] = AllLegs_contactRolling_J_and_Jdot(q, qd, params);


% Initialize total damping torque [Combined rolling + lateral]
tau_damping = DampingInfluence(q,qd,J_Total,params);



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

% Log into parameter object
p.lambda(:, end+1) = lambda;
p.Sys_Input(:, end+1) = Sys_Input;
p.t(end+1) = t;


% Compose derivative of state vector
dy = [qd; qdd];








% === Display simulation progress ===
if toc(lastTime) > 2
    fprintf('Time: %.2f s | Elapsed Time: %.2f s\n', t, toc(lastTime));
    lastTime = tic;
end

end


