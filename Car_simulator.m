clc
clear
close all

%%  Define Vehicle Parameters ===
mb = 150;           % Mass of car body (kg)
mw = 25;            % Mass of each wheel (kg)
Ibx = 400; Iby = 1300; Ibz = 1800;  % Moments of inertia (kg*m^2)
Iw = 2;             % Wheel rotational inertia (kg*m^2)
g = 9.81;           % Gravity (m/s^2)
R = 0.3;            % Wheel radius (m)
l = 1;              % Half length of vehicle (front-rear)
w = 0.5;           % Half width of vehicle (left-right)
lat_damp = -800*1;  % N·s/m, per wheel lateral damping coefficient
roll_damp = -10*1;  % damp_roll (N·s/m) per wheel  [5 for general, 10 for roungh]
n_wheel = 6;

% Pack parameters into a vector matching matlabFunction order
params  = [mb; mw; Ibx; Iby; Ibz; Iw; g; R; l; w; lat_damp; roll_damp; n_wheel];

%%  Initial conditions ===
q0  = zeros(12,1);         % All positions and angles = 0
q0(3) = 0.30; % Height of the chassis
q0(5) = 0.00; % Angle of the chassis
qd0 = zeros(12,1);         % All velocities = 0
qd0(1) = 0.0;

yo = [q0; qd0];  % Combine initial positions and velocities





% ODE function (carDynamics must compute dy = [qd; qdd])
f = @(t, y) carDynamics(t, y, tau_func, params);

% Time simulation setup
tspan = [0 5];



% Get Jacobians
[J,~,~,key_idx] = AllLegs_contactRolling_J_and_Jdot(q0, qd0, params);

[ml,~] = size(J);




% This stores the parameters from inside loop of the ODE/RK 
p = ParamLog;

p.lambda = zeros(ml,1);
p.Sys_Input = zeros(length(q0),1);
p.R = R;
p.w = w;
p.n_wheel = n_wheel;
p.track = 1; % if 0, no tracking is done

%%  Torque Input to the Model ========
% PD controller Torque

% Apply forward torque to all wheels
% tau = @(t, y) 5*[0;0;0;0;0;0;1;1;1;1;1;1];

if p.track == 1
    tau_func = @(t,y) PD_Controller(t, y, p);
else
    % Torque function (adjusted to 6 wheels)
    tau_func = @(t,y) 30*[0;0;0;0;0;0; 1; 1; 1; 1; 1; 1];
end
%%  Simulation Loop  
solver_mode = 'ODE';
% solver_mode = 'RK';

initial_conditions = yo;



%% Parameters for ODE/RK Innerloop setup

% Dynamics Simulations


% reset time [Persistent Variable ]
carDynamics([], [], [], [], [], true);




fprintf("\n Simulating the Model, Please wait .............\n")


switch solver_mode
    case 'ODE'
        fprintf('\n Solving using MATALB ODE SOLVER ........\n');

        F = ode;
        F.Parameters = p;
        F.ODEFcn = @(t, y, p) carDynamics(t, y, p, tau_func, params);

        F.Solver = "ode15s";
        F.SolverOptions.MaxStep = 1e-3; % This may slow down the solution speed
        % F.EventDefinition = E;

        tic
        % profile on;

        F.InitialValue = initial_conditions;
        sol = solve(F,0,tspan(end));
        % profile viewer;
        tout = sol.Time;
        yout = sol.Solution';
        toc;


    case 'RK'
        % RK 4th Order method
        % Using RK Method to get the output
        fprintf('\n Solving using RK 4th order approximation\n........')
        tic
        dt = 0.005;
        tspan_s = 0:dt:tspan(end);
        y0 = initial_conditions;


        yout = zeros(length(y0),length(tspan_s));
        yout(:,1) = y0;

        for i = 1:length(tspan_s)-1

            y = yo;


            % Updated appraoch based on RK method

            t = tspan_s(i);
            % Get lambda at this step using current state
            k1 = carDynamics(t, y, p, tau_func, params);
            k2 = carDynamics(t + dt/2, y + dt/2 * k1, p, tau_func, params);
            k3 = carDynamics(t + dt/2, y + dt/2 * k2, p, tau_func, params);
            k4 = carDynamics(t + dt, y + dt * k3, p, tau_func, params);

            Y1 = y + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);


            % Extract generalized coordinates and velocities from state vector
            q  = y(1:length(q0));
            qd = y(length(q0)+1:20);

            % ===  State Update ====
            yout(:,i+1) = Y1;
            yo = Y1;  % Set initial condition for next step


        end
        yout = yout';
        y = yout';
        tout = tspan_s';
        toc
end

fprintf("\n Model Simulated Successfully. \n\n")








%% Extract Lambda and Forces at Reaction


t_vec = (tout);  % or use t directly if it's a vector
lambda_out = zeros(ml,length(t_vec));  % preallocate
Sys_Input = zeros(length(q0), length(t_vec));  % same columns as t_vec
Tau_in = zeros(p.n_wheel, length(t_vec));  % same columns as t_vec

for i = 1:length(t_vec)
    ti = t_vec(i);
    idx = find(abs(p.t - ti) < 1e-8);

    if ~isempty(idx)

        lambda_out(:,i) = p.lambda(:,idx(min(1,length(idx))));
        Sys_Input(:, i) = p.Sys_Input(:, idx(1));
        if p.track == 1
            Tau_in(:, i) = p.tau(:, idx(1));
        else
            tor = tau_func([],[]);
            Tau_in(:, i) = tor(7:7+p.n_wheel-1);
        end
    else
        lambda_out(:,i) = NaN;  % or keep empty
        Sys_Input(:, i) = NaN;
        Tau_in(:, i) = NaN;
    end
end






%% Plotting


% Extract states groups from yout
pos      = yout(:, 1:3);      % x, y, z
orient   = yout(:, 4:6);      % phi, theta, psi
wheels   = yout(:, 7:12);     % th1, th2, th3, th4
vel_lin  = yout(:, 13:15);    % xd, yd, zd
vel_ang  = yout(:, 16:18);    % phid, thetad, psid
vel_wheel= yout(:, 19:24);    % th1d, th2d, th3d, th4d


%% Plot groups
plotGroup(tout, pos, {'x', 'y', 'z'}, 'Position (m)');

plotGroup(tout, orient, {'$\phi (roll)$', '$\theta (pitch)$', '$\psi (yaw)$'}, 'Orientation (rad)');
%
plotGroup(tout, wheels, {'$\theta_1$', '$\theta_2$', '$\theta_3$', '$\theta_4$', '$\theta_5$', '$\theta_6$'}, 'Wheel Angles (rad)');
%
plotGroup(tout, vel_lin, {'$\dot{x}$', '$\dot{y}$', '$\dot{z}$'}, 'Linear Velocities (m/s)');
%
plotGroup(tout, vel_ang, {'$\dot{\phi}$', '$\dot{\theta}$', '$\dot{\psi}$'}, 'Angular Velocities (rad/s)');
%
plotGroup(tout, vel_wheel, {'$\dot{\theta}_1$', '$\dot{\theta}_2$', '$\dot{\theta}_3$', '$\dot{\theta}_4$', '$\dot{\theta}_5$', '$\dot{\theta}_6$'}, 'Wheel Angular Velocities (rad/s)');



%% Reaction force plotting goes here


figure;

n_f = ml/p.n_wheel;


Force_legends = {'X-direction Contact Forces';
    'Y-direction Contact Forces';
    'Z-direction Contact Forces'};



legendLabels = arrayfun(@(i) sprintf('F_{%d}', i), 1:p.n_wheel, 'UniformOutput', false);



% Subplot 1: X-direction contact forces (odd indices)
for j = 1:n_f
    nexttile
    plot(tout(1:end-1), lambda_out(0+j:n_f:ml,1:end-1).','LineWidth',2);
    xlabel('Time (s)');
    ylabel('F_x (N)');
    title(Force_legends{key_idx(j)});
    legend(legendLabels{:});
    grid on;
end



figure
for j = 1:10
    nexttile()
    plot(Sys_Input(j,1:end-1));
    title(num2str(j));
    axis tight;
    xlabel('Time (s)');
end




%% Plot the torques at the wheel

legendLabels = arrayfun(@(i) sprintf('T_{%d}', i), 1:p.n_wheel, 'UniformOutput', false);

figure;

plot(tout(1:end-1), Tau_in(:,1:end-1).','LineWidth',2);
xlabel('Time (s)');
ylabel('Torque (Nm)');
title('Torque Input at Wheels');
legend(legendLabels{:});
grid on;



%% Animation of motion



% VehicleMotion_Animation


