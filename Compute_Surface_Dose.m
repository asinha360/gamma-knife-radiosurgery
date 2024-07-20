function [points, point_doses] = Compute_Surface_Dose(num_points, sphere_center, sphere_radius, beam_nvectors, isocenter, skin_entry_points, dose_absorption_function_table, radial_dose_function_table, indices_unsafe)
    % Compute_Surface_Dose generates points on the surface of a sphere and computes the dose at each point.
    % Inputs:
    %   num_points - Number of points to generate on the sphere's surface
    %   sphere_center - Center of the sphere
    %   sphere_radius - Radius of the sphere
    %   beam_nvectors - Array of beam direction vectors
    %   isocenter - The isocenter of the beams
    %   skin_entry_points - Array of skin entry points for each beam
    %   dose_absorption_function_table - Table containing depth dose information
    %   radial_dose_function_table - Table containing radial dose information
    %   indices_unsafe - Indices of beams considered unsafe
    % Outputs:
    %   points - Array of points on the sphere's surface
    %   point_doses - Array of calculated doses at each point

    % Initialize the seed for random number generation
    rng(0, 'twister');

    % Generate random points on the surface of the sphere
    rvals = 2 * rand(num_points, 1) - 1;
    elevation = asin(rvals);
    azimuth = 2 * pi * rand(num_points, 1);
    radii = sphere_radius * ones(num_points, 1);
    [x, y, z] = sph2cart(azimuth, elevation, radii);

    % Offset the points by the sphere's center
    X = x + sphere_center(1);
    Y = y + sphere_center(2);
    Z = z + sphere_center(3);

    % Initialize arrays for storing points and their corresponding doses
    points = zeros(num_points, 3);
    point_doses = zeros(num_points, 1);

    % Compute dose for each point on the sphere's surface
    for i = 1:num_points
        point_of_interest = [X(i), Y(i), Z(i)];
        points(i, :) = point_of_interest;

        % Calculate the dose at each point considering all beams
        point_doses(i) = Compute_Dose_from_All_Beams(point_of_interest, beam_nvectors, isocenter, skin_entry_points, dose_absorption_function_table, radial_dose_function_table, indices_unsafe);
    end
end