clc;
clear;
close all;

% === Define symbols ===
syms x y z phi theta psi real   % COM position and orientation
syms xd yd zd phid thetad psid real
syms xdd ydd zdd phidd thetadd psidd real
syms th1 th2 th3 th4 th5 th6 real
syms th1d th2d th3d th4d th5d th6d real
syms th1dd th2dd th3dd th4dd th5dd th6dd real
syms mb mw Ibx Iby Ibz Iw g R l w real
syms tau1 tau2 tau3 tau4 tau5 tau6 real

% Generalized coordinates and derivatives
q = [x; y; z; phi; theta; psi; th1; th2; th3; th4; th5; th6];
qd = [xd; yd; zd; phid; thetad; psid; th1d; th2d; th3d; th4d; th5d; th6d];
qdd = [xdd; ydd; zdd; phidd; thetadd; psidd; th1dd; th2dd; th3dd; th4dd; th5dd; th6dd];
nq = length(q);

vars_q  = q;
vars_qd = qd;
params  = [mb; mw; Ibx; Iby; Ibz; Iw; g; R; l; w];

% Transformation from Euler angle rates to body angular velocity
T_euler = [1, 0, -sin(theta);
           0, cos(phi), sin(phi)*cos(theta);
           0, -sin(phi), cos(phi)*cos(theta)];
omega_body = simplify(T_euler * [phid; thetad; psid]);

% Inertia matrix of the car body
Ib = diag([Ibx, Iby, Ibz]);

% Kinetic energy of the body
T_body = (1/2)*mb*(xd^2 + yd^2 + zd^2) + (1/2)*(omega_body.' * Ib * omega_body);

% Wheel hub positions
r_FL = [ l;  w; 0];
r_FR = [ l; -w; 0];
r_ML = [ 0;  w; 0];
r_MR = [ 0; -w; 0];
r_RL = [-l;  w; 0];
r_RR = [-l; -w; 0];

r_FL_contact_b = r_FL + [0; 0; -R];
r_FR_contact_b = r_FR + [0; 0; -R];
r_ML_contact_b = r_ML + [0; 0; -R];
r_MR_contact_b = r_MR + [0; 0; -R];
r_RL_contact_b = r_RL + [0; 0; -R];
r_RR_contact_b = r_RR + [0; 0; -R];

p_body = [x; y; z];
R_body = combinedRotationMatrix(phi, theta, psi);

r_FL_contact_w = R_body * r_FL_contact_b + p_body;
r_FR_contact_w = R_body * r_FR_contact_b + p_body;
r_ML_contact_w = R_body * r_ML_contact_b + p_body;
r_MR_contact_w = R_body * r_MR_contact_b + p_body;
r_RL_contact_w = R_body * r_RL_contact_b + p_body;
r_RR_contact_w = R_body * r_RR_contact_b + p_body;

J_FL = jacobian(r_FL_contact_w, q);
J_FR = jacobian(r_FR_contact_w, q);
J_ML = jacobian(r_ML_contact_w, q);
J_MR = jacobian(r_MR_contact_w, q);
J_RL = jacobian(r_RL_contact_w, q);
J_RR = jacobian(r_RR_contact_w, q);

% Wheel angle selector vectors
e_theta1 = zeros(1,nq); e_theta1(7) = 1;
e_theta2 = zeros(1,nq); e_theta2(8) = 1;
e_theta3 = zeros(1,nq); e_theta3(9) = 1;
e_theta4 = zeros(1,nq); e_theta4(10) = 1;
e_theta5 = zeros(1,nq); e_theta5(11) = 1;
e_theta6 = zeros(1,nq); e_theta6(12) = 1;

e_rolling_body = [1; 0; 0];
e_vertical_body = [0; 0; 1];

e_rolling = R_body * e_rolling_body;
e_vertical = R_body * e_vertical_body;

Jc_FL = [ e_rolling.' * J_FL - R * e_theta1; e_vertical.' * J_FL ];
Jc_FR = [ e_rolling.' * J_FR - R * e_theta2; e_vertical.' * J_FR ];
Jc_ML = [ e_rolling.' * J_ML - R * e_theta5; e_vertical.' * J_ML ];
Jc_MR = [ e_rolling.' * J_MR - R * e_theta6; e_vertical.' * J_MR ];
Jc_RL = [ e_rolling.' * J_RL - R * e_theta3; e_vertical.' * J_RL ];
Jc_RR = [ e_rolling.' * J_RR - R * e_theta4; e_vertical.' * J_RR ];

Jc = [Jc_FL; Jc_FR; Jc_ML; Jc_MR; Jc_RL; Jc_RR];

%% === Compute time derivatives Jdot_c ===
Jdot_c_FL = sym(zeros(size(Jc_FL)));
Jdot_c_FR = sym(zeros(size(Jc_FR)));
Jdot_c_ML = sym(zeros(size(Jc_ML)));
Jdot_c_MR = sym(zeros(size(Jc_MR)));
Jdot_c_RL = sym(zeros(size(Jc_RL)));
Jdot_c_RR = sym(zeros(size(Jc_RR)));

[m,n] = size(Jc_FL);

for row = 1:m
    for i = 1:n
        dJ_dq1 = diff(Jc_FL(row,:), q(i));
        Jdot_c_FL(row,:) = Jdot_c_FL(row,:) + dJ_dq1 * qd(i);

        dJ_dq2 = diff(Jc_FR(row,:), q(i));
        Jdot_c_FR(row,:) = Jdot_c_FR(row,:) + dJ_dq2 * qd(i);

        dJ_dq3 = diff(Jc_ML(row,:), q(i));
        Jdot_c_ML(row,:) = Jdot_c_ML(row,:) + dJ_dq3 * qd(i);

        dJ_dq4 = diff(Jc_MR(row,:), q(i));
        Jdot_c_MR(row,:) = Jdot_c_MR(row,:) + dJ_dq4 * qd(i);

        dJ_dq5 = diff(Jc_RL(row,:), q(i));
        Jdot_c_RL(row,:) = Jdot_c_RL(row,:) + dJ_dq5 * qd(i);

        dJ_dq6 = diff(Jc_RR(row,:), q(i));
        Jdot_c_RR(row,:) = Jdot_c_RR(row,:) + dJ_dq6 * qd(i);
    end
end



Jdot_c = [Jdot_c_FL;
            Jdot_c_FR;
            Jdot_c_ML;
            Jdot_c_MR;
            Jdot_c_RL;
            Jdot_c_RR];   % 12 x 10
%% === Export MATLAB functions ===


matlabFunction(Jc, Jdot_c, ...
    'Vars', {vars_q, vars_qd, params}, ...
    'File', 'AllLegs_contactRolling_J_and_Jdot', ...
    'Outputs', {'Jc', 'Jdot_c'});

















v_com = [xd; yd; zd];
v_FL = v_com + cross(omega_body, r_FL);
v_FR = v_com + cross(omega_body, r_FR);
v_ML = v_com + cross(omega_body, r_ML);
v_MR = v_com + cross(omega_body, r_MR);
v_RL = v_com + cross(omega_body, r_RL);
v_RR = v_com + cross(omega_body, r_RR);

T_wheel_trans = (1/2)*mw*(v_FL.'*v_FL + v_FR.'*v_FR + v_ML.'*v_ML + v_MR.'*v_MR + v_RL.'*v_RL + v_RR.'*v_RR);
T_wheel_rot   = (1/2)*Iw*(th1d^2 + th2d^2 + th3d^2 + th4d^2 + th5d^2 + th6d^2);

T_total = simplify(T_body + T_wheel_trans + T_wheel_rot);

z_FL = p_body + (R_body * r_FL);
z_FR = p_body + (R_body * r_FR);
z_ML = p_body + (R_body * r_ML);
z_MR = p_body + (R_body * r_MR);
z_RL = p_body + (R_body * r_RL);
z_RR = p_body + (R_body * r_RR);

V_total = mb * g * z + mw * g * (z_FL(3) + z_FR(3) + z_ML(3) + z_MR(3) + z_RL(3) + z_RR(3));

L = T_total - V_total;
Tau = [zeros(6,1); tau1; tau2; tau3; tau4; tau5; tau6];

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

%%
% EXporting the functions

matlabFunction(M, 'Vars', {vars_q, params}, 'File', 'M_matrix_func');
matlabFunction(C, 'Vars', {vars_q, vars_qd, params}, 'File', 'C_vector_func');
matlabFunction(G, 'Vars', {vars_q, params}, 'File', 'G_vector_func');



