function [dose_absorption_function_table] = Compute_Depth_Dose(max_head_size, increment)
    % Compute_Depth_Dose builds a lookup table for the Dose Absorption Function
    % from 0 to the maximum head size.
    % Input:
    %   max_head_size - The maximum size of the head in mm.
    %   increment - The increment in depth in mm.
    % Output:
    %   dose_absorption_function_table - A table with columns 'Depth' in mm
    %   and 'Dose' in %.

    % Calculate the number of depth points
    num_depth_points = floor(max_head_size / increment) + 1;

    % Initialize Depth and Dose arrays
    Depth = (0:increment:max_head_size)';
    Dose = zeros(num_depth_points, 1);

    % Compute dose values based on skin depth
    for i = 1:num_depth_points
        depth_from_skin = Depth(i);
        if depth_from_skin <= 20
            Dose(i) = 0.025 * depth_from_skin + 0.5;
        else
            Dose(i) = -0.005 * depth_from_skin + 1.1;
        end
    end

    % Build table
    dose_absorption_function_table = table(Depth, Dose, 'VariableNames', {'Depth', 'Dose'});
end