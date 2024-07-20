function [points, point_doses] = Compute_Dose(dose_box, increment, beam_nvectors, isocenter, skin_entry_points, dose_absorption_function_table, radial_dose_function_table, indices_unsafe)
    % Compute_Dose generates equidistant points within a specified dose box and computes the dose at each point.
    % Inputs:
    %   dose_box - A 1x6 array defining the lower and upper corners of the dose box [x_min, y_min, z_min, x_max, y_max, z_max]
    %   increment - The spacing between points in the dose box
    %   beam_nvectors - Array of beam direction vectors
    %   isocenter - Central point of the radiation beams
    %   skin_entry_points - Entry points of beams on the skin surface
    %   dose_absorption_function_table - Table for depth dose calculations
    %   radial_dose_function_table - Table for radial dose calculations
    %   indices_unsafe - Indices of beams considered unsafe
    % Outputs:
    %   points - Array of generated points within the dose box
    %   point_doses - Corresponding dose values at each point

    % Determine the number of grid points in x, y, and z directions
    numx = length(dose_box(1)-dose_box(4)/2 : increment : dose_box(1)+dose_box(4)/2);
    numy = length(dose_box(2)-dose_box(5)/2 : increment : dose_box(2)+dose_box(5)/2);
    numz = length(dose_box(3)-dose_box(6)/2 : increment : dose_box(3)+dose_box(6)/2);

    % Initialize arrays for storing points and doses
    total_points = numx * numy * numz;
    points = zeros(total_points, 3);
    point_doses = zeros(total_points, 1);

    % Generate points and compute dose at each point
    ind = 1; % Index for storing results
    for xbox = dose_box(1)-dose_box(4)/2 : increment : dose_box(1)+dose_box(4)/2
        for ybox = dose_box(2)-dose_box(5)/2 : increment : dose_box(2)+dose_box(5)/2
            for zbox = dose_box(3)-dose_box(6)/2 : increment : dose_box(3)+dose_box(6)/2
                point_of_interest = [xbox, ybox, zbox];
                points(ind, :) = point_of_interest;

                % Compute the dose at the current point
                point_doses(ind) = Compute_Dose_from_All_Beams(point_of_interest, beam_nvectors, isocenter, skin_entry_points, dose_absorption_function_table, radial_dose_function_table, indices_unsafe);
                ind = ind + 1;
            end
        end
    end
end