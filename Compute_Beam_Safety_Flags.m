function [beam_safety_flags, num_unsafe_beams] = Compute_Beam_Safety_Flags(oar_center, oar_radius, beam_diameter, isocenter, beam_nvectors)
    % Compute_Beam_Safety_Flags determines the safety flags for each beam.
    % A beam is considered unsafe if it intersects at all with the OAR (Organ At Risk).
    % Input:
    %   oar_center, oar_radius - Center and radius of the OAR
    %   beam_diameter - Diameter of the beam
    %   isocenter - The target isocenter
    %   beam_nvectors - Direction vectors of the beams
    % Output:
    %   beam_safety_flags - A binary 1D array where 1 indicates a beam intersecting with the OAR
    %   num_unsafe_beams - The number of beams that intersect with the OAR

    % Initialize array for intersections
    numIntersections = zeros(size(beam_nvectors, 1), 1);

    % Check intersections for each beam
    for i = 1:length(beam_nvectors)
        numIntersections(i) = Num_Intersections_Of_Sphere_And_Cylinder(oar_center, oar_radius, beam_diameter / 2, isocenter, -beam_nvectors(i, :));
    end

    % Determine safety flags (1 for unsafe beams)
    beam_safety_flags = numIntersections > 0;

    % Compute number of unsafe beams
    num_unsafe_beams = sum(beam_safety_flags);
end