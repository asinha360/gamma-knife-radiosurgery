function [beam_nvectors] = Compute_Beam_Directions(beam_separation_angle)
    % Compute_Beam_Directions computes the unit direction vector for each pencil
    % beam's centerline.
    % Input:
    %   beam_separation_angle - The angle between beams in degrees.
    % Output:
    %   beam_nvectors - A 2D array with each row corresponding to the [x,y,z] unit vector directions.

    % Calculate the number of rotations
    num_rot_x = round(90 / beam_separation_angle);
    num_rot_z = round(360 / beam_separation_angle);

    % Preallocate beam vector array
    total_vectors = num_rot_x * num_rot_z + 1;
    beam_nvectors = zeros(total_vectors, 3);

    % Start with z unit vector
    z_unit_vector = [0, 0, 1];

    % Initialize index for beam vector array
    index = 2;
    beam_nvectors(1, :) = z_unit_vector; % Add first z axis vector

    % Rotate about x and z by angle increments
    for i = 0:num_rot_x - 1
        anglex = (i + 1) * beam_separation_angle; % Degrees
        [~, R4by4x] = Rotation_About_Frame_Axis(1, anglex);

        for j = 0:num_rot_z - 1
            anglez = j * beam_separation_angle; % Degrees
            [~, R4by4z] = Rotation_About_Frame_Axis(3, anglez);

            % Rotate z unit vector by all combinations of rotations in x and z
            nvector = (R4by4z * R4by4x * [z_unit_vector, 1]')';
            beam_nvectors(index, :) = nvector(1:3);
            index = index + 1;
        end
    end
end