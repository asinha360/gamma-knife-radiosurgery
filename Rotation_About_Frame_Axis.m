function [R3by3, R4by4] = Rotation_About_Frame_Axis(axis, angle)
    % Rotation_About_Frame_Axis generates a transformation matrix to rotate
    % a point about one of the (x, y, z) frame axes by a given rotation angle
    % in both 3x3 and 4x4 form.
    % Input: 
    %   axis - 1, 2, or 3 for x, y, z respectively
    %   angle - rotation angle in degrees, counterclockwise when looking in the negative direction
    % Output: 
    %   R3by3 - 3x3 rotation matrix
    %   R4by4 - 4x4 homogeneous rotation matrix

    % Convert angle from degrees to radians
    t = deg2rad(angle);

    % Initialize 4x4 homogeneous rotation matrix
    R4by4 = eye(4);

    % Compute rotation matrices based on the specified axis
    switch axis
        case 1 % Rotation about x-axis
            R3by3 = [1, 0, 0; 0, cos(t), -sin(t); 0, sin(t), cos(t)];
        case 2 % Rotation about y-axis
            R3by3 = [cos(t), 0, sin(t); 0, 1, 0; -sin(t), 0, cos(t)];
        case 3 % Rotation about z-axis
            R3by3 = [cos(t), -sin(t), 0; sin(t), cos(t), 0; 0, 0, 1];
        otherwise
            error('Invalid axis. Axis must be 1, 2, or 3 (for x, y, z).');
    end

    % Populate the upper-left 3x3 part of the 4x4 matrix with the 3x3 rotation matrix
    R4by4(1:3, 1:3) = R3by3;
end