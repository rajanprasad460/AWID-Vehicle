function tau_out = PD_Controller(t, y, p)
% PD_Controller - Computes wheel torques using a PD control strategy
% and a torque-to-force allocation system.
%
% Inputs:
%   ~  - Placeholder for time (t), unused here.
%   y  - Current system state vector [q; qd].
%   p  - Parameter structure containing:
%           R      - Wheel radius
%           w      - Half track width
%           tau0   - Initial torque guess for optimizer
%
% Output:
%   tau_out - Full torque vector [0 0 0 0 0 0 tau_wheels]






% === Extract parameters from p ===
R       = p.R;         % Wheel radius
w       = p.w;         % Half track width
tau0    = p.tau0;      % Initial guess for torque distribution
n_wheel = p.n_wheel;           % Number of wheels

% === State unpacking ===
q  = y(1:length(y)/2);     % Generalized coordinates
qd = y(length(y)/2+1:end); % Generalized velocities

% === Wheel kinematics: torque-to-force mapping ===
% Converts longitudinal forces & yaw moment into wheel torques
% A maps [Fx; Mz] to wheel torques: tau_wheels = A \ [Fx; Mz]
A = (1 / (R)) * ...
    [ 1,   1,  1,  1, 1, 1;    % Longitudinal force contribution
    w,  -w,  w, -w , w, -w];  % Yaw moment contribution

% Torque limits (N·m)
lb = -215 * ones(n_wheel, 1);
ub =  215 * ones(n_wheel, 1);

% === Body rotation ===
% Transform global velocities into body frame
R_body = combinedRotationMatrix(q(4), q(5), q(6));  % Roll, pitch, yaw
% R_body = eye(3);
v_body = R_body * qd(1:3);    % Linear velocity in body frame
w_body = R_body * qd(4:6);    % Angular velocity in body frame


% === Control activation switches ===
t_vx     = 0;  % 1 = enable forward/reverse velocity control
t_yaw_z  = 1;  % 1 = enable yaw rate control

% === Velocity tracking errors (body-x) ===
ex   = (1.0 - v_body(1)) * t_vx;     % velocity error (m/s)
vx_accel_est = (v_body(1) - p.dy(13))*0;
edx  = (0.0 - vx_accel_est) * t_vx;  % accel error (m/s²)

% === Yaw rate tracking errors (body-z) ===
e_psi = (-0.2 + w_body(3)) * t_yaw_z; % yaw rate error (rad/s)
yaw_accel_est = (w_body(3) - p.dy(18))*0;
ed_psi = (0.0 + yaw_accel_est) * t_yaw_z; % yaw accel error (rad/s²)





% === PD Control Gains (tuned separately for vx and yaw) ===
Kp = diag([2000, 20000]);   % [Kp_vx, Kp_yaw]
Kd = diag([800, 500]);   % [Kd_vx, Kd_yaw]

% === Combined control output: B = [Fx; Mz] ===
B  = Kp * [ex; e_psi] + Kd * [edx; ed_psi];

% === Torque allocation ===
tau = TractionPlannar(3, A, B, lb, ub, tau0);


% === Store for next call (warm start) ===
p.tau0 = tau';

% === Assemble full torque output vector ===
% First 6 elements are placeholders for chassis/other joints
tau_out = [0 0 0 0 0 0 tau.']';
p.tau(:,end+1) = tau; % Passing tau for visulation after simlulation
end
