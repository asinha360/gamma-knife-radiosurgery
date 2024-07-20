function [doseColors] = Create_Color_Map(x, dose)
    % Create_Color_Map generates a color map based on the provided dose values.
    % This function maps each dose value to a color in a colormap.
    % Inputs:
    %   x - Array representing the positions or indices of dose values.
    %   dose - Array of dose values for which the color map is to be created.
    % Output:
    %   doseColors - Array of colors corresponding to each dose value.

    % Number of points or positions
    numPoints = length(x);

    % Create a colormap using the 'jet' color scheme
    cmap = jet(numPoints);

    % Determine the minimum and maximum dose values
    minDose = min(dose);
    maxDose = max(dose);

    % Calculate the relative position of each dose value within the min-max range
    percentage = (maxDose - dose) / (maxDose - minDose);

    % Determine the index of each dose value in the colormap
    indexes = round(percentage * (numPoints - 1) + 1);

    % Map each dose value to a corresponding color in the colormap
    doseColors = cmap(numPoints - indexes + 1, :);
end