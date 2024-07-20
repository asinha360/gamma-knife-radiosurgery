function [relative_dose, ratio_of_total_structure_volume] = Create_Dose_Volume_Histogram(point_doses, D100, dose_max, inc)
    % Create_Dose_Volume_Histogram calculates the dose volume histogram data.
    % It computes the percentage of the total structure volume that receives at least a certain dose level.
    % Inputs:
    %   point_doses - Array of dose values at each voxel
    %   D100 - Dose value considered as 100% (typically the prescribed dose)
    %   dose_max - Maximum dose value to consider in the histogram
    %   inc - Incremental dose step for histogram calculation
    % Outputs:
    %   relative_dose - Array of dose levels as a percentage of D100
    %   ratio_of_total_structure_volume - Array of percentages of the total volume exceeding each dose level

    % Total number of voxels in the structure
    num_voxels = length(point_doses);

    % Pre-allocate arrays for efficiency
    num_steps = floor(dose_max / inc) + 1;
    ratio_of_total_structure_volume = zeros(num_steps, 1);
    relative_dose = zeros(num_steps, 1);

    % Calculate ratio of total structure volume and relative dose for each dose level
    for idx = 1:num_steps
        current_dose = (idx - 1) * inc;
        ratio_of_total_structure_volume(idx) = sum(point_doses > current_dose) / num_voxels * 100; % Percentage of volume above current dose level
        relative_dose(idx) = current_dose / D100 * 100; % Relative dose as a percentage of D100
    end
end