function depth_from_skin = Compute_Depth_from_Skin(point_of_interest, beam_nvectors, beam_index, isocenter, skin_entry_points)
    % Compute_Depth_from_Skin calculates the depth from the skin to a point of interest along a specific beam.
    % Input:
    %   point_of_interest - The point for which depth is to be calculated
    %   beam_nvectors - Array of beam direction vectors
    %   beam_index - Index of the specific beam to be considered
    %   isocenter - The isocenter point of the beams
    %   skin_entry_points - Array of points where each beam enters the skin
    % Output:
    %   depth_from_skin - The calculated depth from the skin to the point of interest along the specified beam

    % Calculate the radial distance vector from the point of interest to the beam
    [~, radial_distance_vector] = Distance_Between_Line_And_Point(isocenter, isocenter + beam_nvectors(beam_index, :), point_of_interest);

    % Calculate the point on the beam closest to the point of interest
    point_on_beam = point_of_interest + radial_distance_vector;

    % Calculate the depth from the skin entry point to this point on the beam
    depth_from_skin = norm(skin_entry_points(beam_index, :) - point_on_beam);
end