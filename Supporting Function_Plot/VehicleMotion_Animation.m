% === BEGIN MODIFIED ANIMATION CODE FOR SIX-WHEEL VEHICLE WITH TRAILS, TORQUE, EXPORT ===

nw = p.n_wheel;  % Total wheels
wheel_radius = p.R;
wheel_width = 0.15;
car_height = 0.05;
sim_speed = 4;

% save_gif = true;  % Toggle to save as GIF
% save_vid = true; % Toggle to save as mpeg

% save_gif = false;  % Toggle to save as GIF
% save_vid = false; % Toggle to save as mpeg



filename = 'six_wheel_animation';
f_count = 1; % Counting the frame



% Extract trajectory
x_traj = yout(:,1); y_traj = yout(:,2); z_traj = yout(:,3);
pad = 0.8;
x_lim = [min(x_traj)-pad*4*l, max(x_traj)+pad*4*l];
y_lim = [min(y_traj)-pad*4*w, max(y_traj)+pad*4*w];
z_lim = [0 2];

Nwheel = 36;  % Wheel rendering resolution
[wheel_circ_x, wheel_circ_z] = cylinder(wheel_radius, Nwheel);
wheel_circ_y = (ones(size(wheel_circ_x)) .* [0; wheel_width]) - wheel_width/2;

fig = figure('Name', '6-Wheel Vehicle Simulation','units','normalized','outerposition',[0 0 1 1]);
hold on;

% === Store trails ===
wheel_trails = cell(1, nw);
for i = 1:nw
    wheel_trails{i} = [];
end


sim_idx = 1:sim_speed:length(tout);
f_count_max = length(sim_idx);  % <-- Set this to known number of frames

im(f_count_max) = struct('cdata', [], 'colormap', []);  % Preallocate movie frames
im2 = cell(1, f_count_max);                             % Preallocate image cell array


for k_in = 1:length(sim_idx)
    clf;
    hold on;
    k = sim_idx(k_in);

    % Extract states
    x = yout(k,1); y = yout(k,2); z = yout(k,3);
    phi = yout(k,4); theta = yout(k,5); psi = yout(k,6);
    th = yout(k,7:12);  % Wheel angles

    % Rotation matrices
    Rx = [1 0 0; 0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)];
    Ry = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
    Rz = [cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];
    R_body_to_world = Rz * Ry * Rx;

    % Car body top face
    body_corners = [ l,  w, 0;
        l, -w, 0;
        -l, -w, 0;
        -l,  w, 0]';
    body_corners(3,:) = car_height;
    world_body_corners = R_body_to_world * body_corners + [x; y; z];

    % === Wheel centers (FL, FR, ML, MR, RL, RR) ===
    wheel_centers_body = [ l,  w, 0;
        l, -w, 0;
        0,  w, 0;
        0, -w, 0;
        -l,  w, 0;
        -l, -w, 0]';

    for i = 1:nw
        center_world = R_body_to_world * wheel_centers_body(:,i) + [x; y; z];
        wheel_trails{i} = [wheel_trails{i}, center_world];

        [Xc, Zc] = cylinder(wheel_radius, Nwheel);
        Yc = (ones(size(Xc)) .* [0; wheel_width]) - wheel_width/2;
        wheel_local = [Xc(:)'; Yc(:)'; Zc(:)'];

        theta_i = th(i);
        Ry_spin = [cos(theta_i), 0, sin(theta_i);
            0, 1, 0;
            -sin(theta_i), 0, cos(theta_i)];
        wheel_spun = Ry_spin * wheel_local;
        wheel_world = Rz * wheel_spun;

        wheel_world(1,:) = wheel_world(1,:) + center_world(1);
        wheel_world(2,:) = wheel_world(2,:) + center_world(2);
        wheel_world(3,:) = wheel_world(3,:) + center_world(3);

        Xw = reshape(wheel_world(1,:), size(Xc));
        Yw = reshape(wheel_world(2,:), size(Yc));
        Zw = reshape(wheel_world(3,:), size(Zc));

        surf(Xw, Yw, Zw, 'FaceColor', 'k', 'EdgeColor', 'none');

        % Spokes
        spoke_tip_local = [0; 0; -wheel_radius];
        phi_spoke = 2*pi/nw;
        for in = 0:nw-1
            Ry_spin_n = [cos(phi_spoke*in), 0, sin(phi_spoke*in);
                0, 1, 0;
                -sin(phi_spoke*in), 0, cos(phi_spoke*in)];
            tip = Rz * (Ry_spin * Ry_spin_n * spoke_tip_local);
            line = [center_world, center_world + tip];
            spoke_color = 'g';
            plot3(line(1,:), line(2,:), line(3,:), spoke_color, 'LineWidth', 1.5);
        end

        % === Plot wheel trails ===
        plot3(wheel_trails{i}(1,:), wheel_trails{i}(2,:), wheel_trails{i}(3,:), ...
            'Color', [0.4 0.4 0.4], 'LineWidth', 1);
    end

    % === Draw chassis ===
    patch('XData', world_body_corners(1,:), ...
        'YData', world_body_corners(2,:), ...
        'ZData', world_body_corners(3,:), ...
        'FaceColor', 'magenta', 'FaceAlpha', 0.7);

    % === Local axes & torque vector ===
    origin = [x; y; z + car_height];
    L_axis = 0.5;
    quiver3(origin(1), origin(2), origin(3), ...
        R_body_to_world(1,1)*L_axis, R_body_to_world(2,1)*L_axis, R_body_to_world(3,1)*L_axis, ...
        'r', 'LineWidth', 2);
    quiver3(origin(1), origin(2), origin(3), ...
        R_body_to_world(1,2)*L_axis, R_body_to_world(2,2)*L_axis, R_body_to_world(3,2)*L_axis, ...
        'g', 'LineWidth', 2);
    quiver3(origin(1), origin(2), origin(3), ...
        R_body_to_world(1,3)*L_axis, R_body_to_world(2,3)*L_axis, R_body_to_world(3,3)*L_axis, ...
        'b', 'LineWidth', 2);
    plot3(origin(1), origin(2), origin(3), 'ko', 'MarkerFaceColor', 'k');

    % === Torque vector ===
    torque_vector = cross([1;0;0], R_body_to_world(:,1));  % Example: X-axis deviation
    quiver3(origin(1), origin(2), origin(3), ...
        torque_vector(1), torque_vector(2), torque_vector(3), ...
        'm', 'LineWidth', 2, 'MaxHeadSize', 2);

    % Draw ground
    % ground_size = 50;
    % fill3([-ground_size ground_size ground_size -ground_size], ...
    %       [-ground_size -ground_size ground_size ground_size], ...
    %       [0 0 0 0], [0.9 0.9 0.9]);

    % Formatting
    axis equal;
    xlim(x_lim); ylim(y_lim); zlim(z_lim);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(sprintf('6-Wheel Vehicle at t = %.2f s', tout(k)));
    view(30, 20); grid on; lighting gouraud; camlight headlight;
    drawnow;

    % ---------- Save Frame ----------
    frame = getframe(gcf);
    im(f_count) = frame;
    im2{f_count} = frame2im(frame);
    f_count = f_count+ 1;

end


%% ============== Export to GIF /Video ===
% ------------------------- [7] Video Export -------------------------
if save_vid
    f_filename = AutoRename(fullfile(cd,'Output Results'),strcat(filename,'.avi'));

    v = VideoWriter(f_filename,'Motion JPEG AVI');
    v.Quality=100;
    v.FrameRate = 3; % No. of frames per second
    open(v);
    writeVideo(v,im);
    close(v);


    pause(1)
    % Open the created video with player to show frame by frame

    implay(f_filename);  % or 'my_animation.avi'

end

%% ------------------------- [9] GIF Export -------------------------

if save_gif

    o_f_name = AutoRename(fullfile(cd,'Output Results'),strcat(filename,'.gif'));
    delay = 0.05;                   % Delay time between frames (in seconds)
    for idx = 1:f_count-1
        [A_map, map] = rgb2ind(im2{idx}, 256); % Convert RGB to indexed image

        if idx == 1
            % Create the GIF file
            imwrite(A_map, map, o_f_name, 'gif', 'LoopCount', Inf, 'DelayTime', delay);
        else
            % Append to the existing GIF
            imwrite(A_map, map, o_f_name, 'gif', 'WriteMode', 'append', 'DelayTime', delay);
        end
    end

end

