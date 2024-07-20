function [skin_entry_points, depth_to_isocenter] = Compute_Skin_Entry_Points(beam_nvectors, isocenter, head_a, head_b, head_c)
    % Compute_Skin_Entry_Points Computes the skin entry point for each beam as
    % well as the depth that the isocenter is from that skin entry point.
    % Input:
    %   beam_nvectors - array of beam direction vectors
    %   isocenter - a 3D point at the center of the PTV
    %   head_a, head_b, head_c - half lengths of the principal axes of an ellipsoid
    % Output:
    %   skin_entry_points - 2D array where each row is the [x,y,z] skin entry point
    %   depth_to_isocenter - 1D array of Euclidean distances from each skin entry point to the isocenter
    
    %Initialize arrays
    skin_entry_points = zeros(size(beam_nvectors));
    depth_to_isocenter = zeros(length(beam_nvectors),1);
    
    for i = 1:length(beam_nvectors)
        %Compute intersection point of beam with head
        [~,skinpoint,~] = Intersection_Of_Vector_With_Ellipsoid(isocenter,-beam_nvectors(i,:),head_a,head_b,head_c);
        skin_entry_points(i,:) = skinpoint;   
        
        %Compute distance between skin entry point and isocenter
        depth_to_isocenter(i) = pdist([skin_entry_points(i,:);isocenter],'euclidean'); 
    end
end