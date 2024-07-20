function [dose_box] = Compute_Dose_Box(ptv_center, ptv_radius, oar_center, oar_radius)
    % Compute_Dose_Box calculates the bounding box encompassing both the 
    % Prescribed Target Volume (PTV) and the Organ At Risk (OAR).
    % Input:
    %   ptv_center, oar_center - 1x3 arrays representing the 3D center coordinates (in mm)
    %   ptv_radius, oar_radius - Scalars representing the radii of PTV and OAR (in mm)
    % Output:
    %   dose_box - 1x6 array containing the 3D coordinates of the bounding box center 
    %   and its side lengths (in mm)

    % Calculate the farthest extents in each direction
    ptv_extents = [ptv_center + ptv_radius; ptv_center - ptv_radius];
    oar_extents = [oar_center + oar_radius; oar_center - oar_radius];
    
    % Combine and find the overall extents
    combined_extents = [ptv_extents; oar_extents];
    box_max = max(combined_extents);
    box_min = min(combined_extents);

    % Calculate box center and side lengths
    box_center = (box_max + box_min) / 2;
    box_side_lengths = box_max - box_min;

    % Output dose box
    dose_box = [box_center, box_side_lengths];
end