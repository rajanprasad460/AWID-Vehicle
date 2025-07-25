clc
clear
close all

% === Define parameters ===
mb = 150;           % Mass of car body (kg)
mw = 25;            % Mass of each wheel (kg)
Ibx = 400; Iby = 1300; Ibz = 1800;  % Moments of inertia (kg*m^2)
Iw = 2;             % Wheel rotational inertia (kg*m^2)
g = 9.81;           % Gravity (m/s^2)
R = 0.3;            % Wheel radius (m)
l = 1.2;            % Half length of vehicle (front-rear)
w = 0.75;           % Half width of vehicle (left-right)
lat_damp = -200;  % N·s/m, per wheel lateral damping coefficient 
n_wheel = 6;

% ====== Spoke parameters ======
nw = 6;      % Number of spokes
% Geometry of spokes
ang_w = 0:2*pi/nw:(2*pi - 2*pi/nw);  % Angles between spokes (rad)
ang_d = 0:360/nw:(360 - 360/nw);    % Degrees between spokes
phi = 2*pi/nw;

% === Initial conditions ===
q0  = zeros(12,1);         % All positions and angles = 0
q0(3) = 0.30; % Height of the chassis
q0(5) = 0.00; % Angle of the chassis
qd0 = zeros(12,1);         % All velocities = 0
qd0(1) = 0.0;

yo = [q0; qd0];  % Combine initial positions and velocities

% Pack parameters into a vector matching matlabFunction order
params  = [mb; mw; Ibx; Iby; Ibz; Iw; g; R; l; w; lat_damp; n_wheel];

% Torque function (adjusted to 6 wheels)
tau_func = @(t) 10*[0;0;0;0;0;0; 1; -0.1; 1; -1; 1; -0.1];

% ODE function (carDynamics must compute dy = [qd; qdd])
f = @(t, y) carDynamics(t, y, tau_func, params);

% Time simulation setup
dt = 0.01;
t_end = 10;
ts = 0:dt:t_end;
ia = length(ts);
yout = zeros(length(ts),24);

% Angular velocity switching
av_mode = 1;

% Initial state
yout(1,:) = yo;

% Get Jacobians
[J,~,~] = AllLegs_contactRolling_J_and_Jdot(q0, qd0, params);
[ml,~] = size(J);
lambda_out = zeros(ml,ia);
Sys_Input = zeros(ml,ia);

% ------------------------- [3] Simulation Loop -------------------------
for i = 1:ia-1
    y = yo;
    t = ts(i);

    [k1, lambda_out(:,i), Sys_Input(:,i)] = carDynamics(t, y, tau_func, params);
    k2 = carDynamics(t + dt/2, y + dt/2 * k1, tau_func, params);
    k3 = carDynamics(t + dt/2, y + dt/2 * k2, tau_func, params);
    k4 = carDynamics(t + dt, y + dt * k3, tau_func, params);

    Y1 = y + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);

    q  = y(1:12);
    qd = y(13:24);

    % % === Spoke impacts ====
    % for iwc = 1:6
    %     if Y1(6+iwc) >= (phi/2)
    %         switch av_mode
    %             case 1
    %                 Y1(6+iwc) = Y1(6+iwc) - phi;
    %                 % Y1(12+iwc) = Y1(12+iwc) * cos(phi);
    %             case 2
    %                 Y1(6+iwc) = Y1(6+iwc) - phi;
    %                 td = Theta_Plus(q,qd,params);
    %                 Y1(12+iwc) = td(iwc);
    %         end
    %     elseif Y1(6+iwc) <= (-phi/2)
    %         switch av_mode
    %             case 1
    %                 Y1(6+iwc) = Y1(6+iwc) + phi;
    %                 % Y1(12+iwc) = Y1(12+iwc) * cos(phi);
    %             case 2
    %                 Y1(6+iwc) = Y1(6+iwc) + phi;
    %                 td = Theta_Plus(q,qd,params);
    %                 Y1(12+iwc) = td(iwc);
    %         end
    %     end
    % end

    yout(i+1,:) = Y1;
    yo = Y1;
end

tout = ts;

%% Extract and plot
pos      = yout(:, 1:3);
orient   = yout(:, 4:6);
wheels   = yout(:, 7:12);
vel_lin  = yout(:, 13:15);
vel_ang  = yout(:, 16:18);
vel_wheel= yout(:, 19:24);

plotGroup(tout, pos, {'x', 'y', 'z'}, 'Position (m)');
plotGroup(tout, orient, {'$\phi$', '$\theta$', '$\psi$'}, 'Orientation (rad)');
plotGroup(tout, wheels, {'$\theta_1$', '$\theta_2$', '$\theta_3$', '$\theta_4$', '$\theta_5$', '$\theta_6$'}, 'Wheel Angles (rad)');
plotGroup(tout, vel_lin, {'$\dot{x}$', '$\dot{y}$', '$\dot{z}$'}, 'Linear Velocities (m/s)');
plotGroup(tout, vel_ang, {'$\dot{\phi}$', '$\dot{\theta}$', '$\dot{\psi}$'}, 'Angular Velocities (rad/s)');
plotGroup(tout, vel_wheel, {'$\dot{\theta}_1$', '$\dot{\theta}_2$', '$\dot{\theta}_3$', '$\dot{\theta}_4$', '$\dot{\theta}_5$', '$\dot{\theta}_6$'}, 'Wheel Angular Velocities (rad/s)');

%% Reaction force plotting
figure;
n_f = ml/6;
for j = 1:n_f
    nexttile
    plot(tout(1:end-1), lambda_out(0+j:n_f:ml,1:end-1).','LineWidth',2);
    xlabel('Time (s)'); ylabel('F_x (N)'); title('X-direction Contact Forces');
    legend('F_1','F_2','F_3','F_4','F_5','F_6');
    grid on;
end

figure
for j = 1:12
    nexttile()
    plot(Sys_Input(j,1:end-1));
    title(num2str(j));
    axis tight;
    xlabel('Time (s)');
end
