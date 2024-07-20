function full_point_dose_value = Compute_Dose_from_All_Beams(point_of_interest, beam_nvectors, isocenter, skin_entry_points, dose_absorption_function_table, radial_dose_function_table, indices_unsafe)
    % Compute_Dose_from_All_Beams computes the cumulative dose at a given point 
    % of interest from all the beams, excluding those identified as unsafe.
    % Input:
    %   point_of_interest - The point for which the dose is being calculated
    %   beam_nvectors - Array of beam direction vectors
    %   isocenter - The isocenter of the beams
    %   skin_entry_points - Array of skin entry points for each beam
    %   dose_absorption_function_table - Table containing depth dose information
    %   radial_dose_function_table - Table containing radial dose information
    %   indices_unsafe - Indices of beams considered unsafe
    % Output:
    %   full_point_dose_value - The total dose contributed by all safe beams at the point of interest

    % Initialize an array to hold dose values from each beam
    point_dose_values = zeros(length(beam_nvectors), 1);

    % Iterate through each beam to calculate its dose contribution
    for i = 1:length(beam_nvectors)
        if ismember(i, indices_unsafe)
            continue; % Skip dose calculation for unsafe beams
        end

        % Calculate radial distance and depth from skin for the current beam
        radial_distance = Compute_Radial_Distance(point_of_interest, beam_nvectors, i, isocenter);
        depth_from_skin = Compute_Depth_from_Skin(point_of_interest, beam_nvectors, i, isocenter, skin_entry_points);

        % Calculate and store the dose contribution from this beam
        point_dose_values(i) = Compute_Point_Dose_from_Beam(dose_absorption_function_table, depth_from_skin, radial_dose_function_table, radial_distance);
    end

    % Sum up the dose contributions from all safe beams
    full_point_dose_value = sum(point_dose_values);
end