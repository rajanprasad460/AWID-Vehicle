function R = combinedRotationMatrix(phi, theta, psi)
%# codegen
    % Rotation matrix around the x-axis (Roll)
    Rx = [1, 0, 0; ...
          0, cos(phi), -sin(phi); ...
          0, sin(phi), cos(phi)];
    
    % Rotation matrix around the y-axis (Pitch)
    Ry = [cos(theta), 0, sin(theta); ...
          0, 1, 0; ...
          -sin(theta), 0, cos(theta)];
    
    % Rotation matrix around the z-axis (Yaw)
    Rz = [cos(psi), -sin(psi), 0; ...
          sin(psi), cos(psi), 0; ...
          0, 0, 1];
    
    % Combined rotation matrix: R = Rx * Ry * Rz
    R = Rz*Ry*Rx;


end