function [distance, dvector] = Distance_Between_Line_And_Point(L1, L2, P)
    % Distance_Between_Line_And_Point calculates the shortest distance between a line and a point.
    % Inputs:
    %   L1, L2 - Two points defining the line
    %   P - The point for which the shortest distance to the line is being calculated
    % Outputs:
    %   distance - The shortest distance from the point P to the line defined by L1 and L2
    %   dvector - The vector from point P to the closest point on the line

    % Vector from the first point on the line (L1) to the point P
    vector_L1_to_P = P - L1;

    % Direction vector of the line from L1 to L2
    line_direction = L2 - L1;
    normalized_line_direction = line_direction / norm(line_direction);

    % Project vector_L1_to_P onto the line's direction vector
    projection = dot(vector_L1_to_P, normalized_line_direction) * normalized_line_direction;

    % Calculate the vector from P to the nearest point on the line
    dvector = vector_L1_to_P - projection;

    % Calculate the magnitude of the distance vector
    distance = norm(dvector);
end