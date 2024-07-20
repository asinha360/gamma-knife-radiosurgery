function point_dose_value = Compute_Point_Dose_from_Beam(dose_absorption_function_table, depth_from_skin, radial_dose_function_table, radial_distance)
    % Compute_Point_Dose_from_Beam computes the dose at a specific point of interest from a beam.
    % Input:
    %   dose_absorption_function_table - Table containing depth dose information
    %   depth_from_skin - Depth from the skin to the point of interest along the beam
    %   radial_dose_function_table - Table containing radial dose information
    %   radial_distance - Radial distance from the beam's central axis to the point of interest
    % Output:
    %   point_dose_value - The calculated dose value at the point of interest

    % Find the closest depth dose value
    [~, ind_closest_depth] = min(abs(dose_absorption_function_table.Depth - depth_from_skin));
    depth_dose = dose_absorption_function_table.Dose(ind_closest_depth);

    % Find the closest radial dose value
    [~, ind_closest_rad] = min(abs(radial_dose_function_table.Radial_Distance - radial_distance));
    radial_dose = radial_dose_function_table.Dose(ind_closest_rad);

    % Calculate the cumulative dose at the point of interest
    point_dose_value = depth_dose * radial_dose;
end