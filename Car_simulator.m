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
lat_damp = -800*0;    % N·s/m, per wheel lateral damping coefficient
roll_damp = -10*1;      % damp_roll (N·s/m) per wheel  [5 for general, 10 for roungh]
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
params  = [mb; mw; Ibx; Iby; Ibz; Iw; g; R; l; w; lat_damp; roll_damp; n_wheel];

% Torque function (adjusted to 6 wheels)
tau_func = @(t) 10*[0;0;0;0;0;0; 1; 1; 1; 1; 1; 1];

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


% =========  Torque to Force conversion ======
A = (1/R)*[1 1 1 1 1 1;w -w w -w w -w];
lb = -500*ones(nw,1);
ub = 500*ones(nw,1);

% F = [F_long;M_psi];
tau0 = 0*ones(1,nw);



% sim_mode = 'ode';
sim_mode = 'RK';

switch sim_mode
    case 'RK'

        % Get Jacobians
        [J,~,~] = AllLegs_contactRolling_J_and_Jdot(q0, qd0, params);
        [ml,~] = size(J);
        lambda_out = zeros(ml,ia);
        Sys_Input = zeros(ml,ia);

        % ------------------------- [3] Simulation Loop -------------------------
        for i = 1:ia-1
            y = yo;
            t = ts(i);

            % === PD Controller =================
            % Errors
            ex = 0;
            ey = 0;
            edx = 2 - yo(13); % Tracking a constant Velocity
            edy = 0 - yo(18); % Yaw rate to  zero

            % Gains
            Kp = 1000;
            Kd = 1000;

            % Desired inertial frame force
            Fx = Kp * ex + Kd * edx;
            My = Kp * ey + Kd * edy;


            B = [Fx;My];
            %---------------   traction optimized distributor ------
            tau = TractionPlannar(3,A,B,lb,ub,tau0);
            tau0 = tau';
            tau_func = @(t) [0 0 0 0 0 0 tau.'].';
            %----------------------------------------


            [k1, lambda_out(:,i), Sys_Input(:,i)] = carDynamics(t, y, tau_func, params);
            k2 = carDynamics(t + dt/2, y + dt/2 * k1, tau_func, params);
            k3 = carDynamics(t + dt/2, y + dt/2 * k2, tau_func, params);
            k4 = carDynamics(t + dt, y + dt * k3, tau_func, params);

            Y1 = y + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);

            q  = y(1:12);
            qd = y(13:24);



            yout(i+1,:) = Y1;
            yo = Y1;


            progressupdater(i,ia-1,'Simulating..');
        end

        tout = ts;
    case 'ode'

        % Reset before starting
        carDynamics_ODE([], [], [], [], true);  % Reset persistent vars

        % Run the simulation
        [tout, yout] = ode15s(@(t,y) carDynamics_ODE(t, y, tau_func, params), ts, yo);


end

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
    xlabel('Time (s)'); ylabel('F (N)'); title(' Contact Forces');
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



%% TO animate

% save_gif = true;  % Toggle to save as GIF
% save_vid = true; % Toggle to save as mpeg

save_gif = false;  % Toggle to save as GIF
save_vid = false; % Toggle to save as mpeg


% VehicleMotion_Animation;