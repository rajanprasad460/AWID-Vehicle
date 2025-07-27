function tau_damping = DampingInfluence(q, qd, J_Total, params)
% DampingInfluence computes the generalized damping torques caused by 
% lateral and rolling damping at the wheels.
%
% Inputs:
%   q        - Generalized coordinates (not directly used here)
%   qd       - Generalized velocities
%   J_Total  - Full Jacobian matrix (3 rows per wheel: [roll; lat; yaw])
%   params   - Vector containing [.., lat_damp, roll_damp, n_wheel]
%
% Output:
%   tau_damping - Generalized damping torque vector (same size as q)

    % Initialize output torque vector
    tau_damping = zeros(length(q), 1);

    % Extract parameters
    lat_damp  = params(11);  % Lateral damping coefficient (N·s/m)
    roll_damp = params(12);  % Rolling damping coefficient (N·s/m)
    n_wheel   = params(13);  % Number of wheels

    for i = 1:n_wheel
        % --- Lateral Damping ---
        % Extract lateral velocity at wheel i using the 2nd row of the Jacobian block
        v_lat_i = J_Total(3*(i-1)+2, :) * qd;
        F_lat_i = -lat_damp * v_lat_i;  % Damping force opposing lateral motion

        % --- Rolling Damping ---
        % Extract rolling velocity at wheel i using the 1st row of the Jacobian block
        v_roll_i = J_Total(3*(i-1)+1, :) * qd;
        F_roll_i = -roll_damp * v_roll_i;  % Damping force opposing rolling motion

        % Combine both damping forces into torque space using Jacobian transpose
        % Ignore yaw damping (3rd row) by setting it to zero
        F_local = [F_roll_i; F_lat_i; 0];  % Local frame damping force/moment vector

        % Accumulate generalized torque
        tau_damping = tau_damping + J_Total(3*(i-1)+1 : 3*i, :)' * F_local;
    end
end
