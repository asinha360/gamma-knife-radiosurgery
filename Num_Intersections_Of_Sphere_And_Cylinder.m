function numInt = Num_Intersections_Of_Sphere_And_Cylinder(C, R, r, P, v)
    % Num_Intersections_Of_Sphere_And_Cylinder computes the number of intersections
    % between a sphere and an infinite cylinder.
    % Returns 0 if there are no intersections, 1 if they just touch,
    % or 2 if the cylinder intersects the sphere at two different locations.
    % Input:
    %   C - Center of the sphere
    %   R - Radius of the sphere
    %   r - Radius of the cylinder
    %   P - A point on the central axis of the cylinder
    %   v - Direction vector of the cylinder's central axis

    % Normalize direction vector of the cylinder
    v = v / norm(v);

    % Find the closest point on the cylinder axis to the sphere center
    t_closest = dot(C - P, v);
    closestPoint = P + t_closest * v;

    % Calculate the distance from the closest point to the sphere center
    d_closest = norm(C - closestPoint);

    % Compute the number of intersections based on geometric relationships
    if d_closest > R + r
        numInt = 0; % No intersection
    elseif d_closest < R - r
        numInt = 2; % Cylinder passes through the sphere
    else
        numInt = 1; % Cylinder just touches the sphere
    end
end