clc;
clear;
close all;

% === Parameters ===
n = 6;  % number of wheels

% === Define symbolic variables ===
syms x y z phi theta psi real
syms xd yd zd phid thetad psid real
syms xdd ydd zdd phidd thetadd psidd real

q = [x; y; z; phi; theta; psi];
qd = [xd; yd; zd; phid; thetad; psid];
qdd = [xdd; ydd; zdd; phidd; thetadd; psidd];

% Wheel symbolic variables
syms mb mw Ibx Iby Ibz Iw g R l w real
params = [mb; mw; Ibx; Iby; Ibz; Iw; g; R; l; w];

q(7:6+n) = 0;
qd(7:6+n) = 0;
qdd(7:6+n) = 0;

for i = 1:n
    syms(['th', num2str(i)], 'real');
    syms(['th', num2str(i), 'd'], 'real');
    syms(['th', num2str(i), 'dd'], 'real');
    syms(['tau', num2str(i)], 'real');

    q(6+i) = eval(['th', num2str(i)]);
    qd(6+i) = eval(['th', num2str(i), 'd']);
    qdd(6+i) = eval(['th', num2str(i), 'dd']);
end

nq = length(q);
vars_q = q;
vars_qd = qd;

% === Rotation Matrix & Angular Velocity ===
T_euler = [1, 0, -sin(theta);
           0, cos(phi), sin(phi)*cos(theta);
           0, -sin(phi), cos(phi)*cos(theta)];
omega_body = simplify(T_euler * qd(4:6));
Ib = diag([Ibx, Iby, Ibz]);

% === Kinetic Energy ===
T_body = (1/2)*mb*sum(qd(1:3).^2) + (1/2)*omega_body.' * Ib * omega_body;

% === Wheel hub positions ===
% Format: [x; y; z] in body frame (predefine for n=6, adjust as needed)
hub_coords = {
    [ l;  w; 0];
    [ l; -w; 0];
    [ 0;  w; 0];
    [ 0; -w; 0];
    [-l;  w; 0];
    [-l; -w; 0];
};

% Transformation
p_body = q(1:3);
R_body = combinedRotationMatrix(phi, theta, psi);
% R_body = eye(3);

tau = sym(zeros(n,1));
T_wheel_trans = 0;
T_wheel_rot = 0;
z_all = sym(zeros(n,1));



e_rolling = R_body * [1; 0; 0];
e_lateral = R_body * [0; 1; 0];
e_vertical = R_body * [0; 0; 1];

e_theta_x = zeros(1,nq);
e_theta_x(1) = 1;


% Filter by index if needed (e.g., only rolling, or all 3 constraints)
key_idx = [1 3];

% Jaocbians

Jc_all = sym(zeros(length(key_idx)*n,nq));
Jdot_c_all = sym(zeros(length(key_idx)*n,nq));
Jc_Total_i = sym(zeros(3,nq));
Jc_Total = sym(zeros(3*n,nq));

for i = 1:n
    r_i = hub_coords{i};
    r_contact_b = r_i + [0; 0; -R];
    r_contact_w = R_body * r_contact_b + p_body;
    
    J_i = jacobian(r_contact_w, q);
    
    e_theta = zeros(1,nq);
    e_theta(6+i) = 1;
    

    Jc_i = [e_rolling.' * J_i - R * e_theta;
            e_lateral.' * J_i;
            e_vertical.' * J_i];
    
    Jc_Total_i = [e_rolling.' * J_i ;
                e_lateral.' * J_i;
                e_vertical.' * J_i];

    Jc_i = Jc_i(key_idx,:);

    % Stack the Jacobians
    Jc_all((i-1)*length(key_idx)+1:(i*length(key_idx)),:) = Jc_i;

    Jc_Total((i-1)*3+1:(i*3),:) = Jc_Total_i;
    
    % Compute Jdot
    mJ = size(Jc_i,1);
    Jdot_i = sym(zeros(mJ,nq));
    for row = 1:mJ
        for j = 1:nq
            dJ = diff(Jc_i(row,:), q(j));
            Jdot_i(row,:) = Jdot_i(row,:) + dJ * qd(j);
        end
    end
    % Jdot_c_all = [Jdot_c_all; Jdot_i];
    % Stack the Jacobians dots
    Jdot_c_all((i-1)*length(key_idx)+1:(i*length(key_idx)),:) = Jdot_i;

    % Velocities
    v_com = qd(1:3);
    v_i = v_com + cross(omega_body, r_i);
    T_wheel_trans = T_wheel_trans + (1/2)*mw*(v_i.'*v_i);
    
    thd = qd(6+i);
    T_wheel_rot = T_wheel_rot + (1/2)*Iw*thd^2;
    
    z_i = p_body + R_body * r_i;
    z_all(i) =  z_i(3);
    
    tau(i) = eval(['tau', num2str(i)]);
end

% === Energies and Lagrangian ===
T_total = simplify(T_body + T_wheel_trans + T_wheel_rot);
V_total = mb*g*q(3) + mw*g*sum(z_all);
L = T_total - V_total;

Tau_vec = [zeros(6,1); tau];

% === Euler-Lagrange Equations ===
EL_eqns = sym(zeros(nq,1));
for i = 1:nq
    dL_dqdot = diff(L, qd(i));
    ddt_dL_dqdot = 0;
    for j = 1:nq
        ddt_dL_dqdot = ddt_dL_dqdot + diff(dL_dqdot, q(j))*qd(j) + diff(dL_dqdot, qd(j))*qdd(j);
    end
    dL_dq = diff(L, q(i));
    EL_eqns(i) = simplify(ddt_dL_dqdot - dL_dq);
end

% === Mass (M), Coriolis (C), Gravity (G) ===
M = sym(zeros(nq));
for i = 1:nq
    for j = 1:nq
        M(i,j) = diff(EL_eqns(i), qdd(j));
    end
end

G = sym(zeros(nq,1));
for i = 1:nq
    G(i) = diff(V_total, q(i));
end

C = simplify(EL_eqns - M*qdd - G);

% === Export Functions ===
matlabFunction(M, 'Vars', {vars_q, params}, 'File', 'M_matrix_func');
matlabFunction(C, 'Vars', {vars_q, vars_qd, params}, 'File', 'C_vector_func');
matlabFunction(G, 'Vars', {vars_q, params}, 'File', 'G_vector_func');

matlabFunction(Jc_all, Jdot_c_all, Jc_Total, key_idx, ...
    'Vars', {vars_q, vars_qd, params}, ...
    'File', 'AllLegs_contactRolling_J_and_Jdot', ...
    'Outputs', {'Jc', 'Jdot_c', 'Jc_Total', 'key_idx'});
