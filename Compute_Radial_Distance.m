function radial_distance = Compute_Radial_Distance(point_of_interest, beam_nvectors, beam_index, isocenter)
    % Compute_Radial_Distance computes the shortest distance from a point of interest
    % to the centerline of a specified beam.
    % Input:
    %   point_of_interest - Coordinates of the point for which the distance is being calculated
    %   beam_nvectors - Array containing the direction vectors of the beams
    %   beam_index - Index specifying which beam's centerline to consider
    %   isocenter - The central point around which the beams are focused
    % Output:
    %   radial_distance - The shortest distance from the point of interest to the beam's centerline

    % Calculate the end point of the beam based on its direction and isocenter
    beam_endpoint = isocenter + beam_nvectors(beam_index, :);

    % Call the function to calculate distance between a line (beam's path) and the point
    radial_distance = Distance_Between_Line_And_Point(isocenter, beam_endpoint, point_of_interest);
end