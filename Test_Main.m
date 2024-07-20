%% CISC 330 A4 Main Function

% Initialization of relevant constants

DOSE_VOXEL_SIZE = 10; % mm
D0 = 1; % mm
D100 = 20; % Prescribed PTV dose units
DOARMAX = 10; % Max OAR dose units

% Helmet
beam_separation_angle = 30; % degrees
beam_diameter = 30; % mm

% Head
head_a = 80; % mm
head_b = 100; % mm
head_c = 80; % mm
head_center = [0, 0, 0]; % mm

% PTV (Prescribed Target Volume) and OAR (Organ At Risk)
ptv_radius = 15; % mm
ptv_center = [30, 0, 15]; % mm
oar_radius = 15; % mm
oar_center = [0, 30, 45]; % mm

isocenter = ptv_center;

%% Draw 3D Scene

fprintf('<strong>Part 1: Drawing the 3-D Gamma Knife Radiosurgery Scene</strong>\n')

% Calculate the "Dose Box" which encompasses the PTV and OAR 
[dose_box] = Compute_Dose_Box(ptv_center, ptv_radius, oar_center, oar_radius);

% Plot the PTV, OAR, Dosebox, Isocenter and the Head
figure;
hold on;
axis equal;  % Maintains the aspect ratio of 3D plot
grid on;     % Enable grid for better orientation

% Plot PTV
[X, Y, Z] = sphere(20);
h1 = surf(ptv_radius*X + ptv_center(1), ptv_radius*Y + ptv_center(2), ptv_radius*Z + ptv_center(3));
set(h1, 'FaceAlpha', 0.4, 'FaceColor', 'cyan', 'EdgeColor', 'none');  % Cyan color for PTV
shading interp;

% Label PTV
text(ptv_center(1) + ptv_radius, ptv_center(2), ptv_center(3), 'PTV', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 8);

% Plot OAR
h2 = surf(oar_radius*X + oar_center(1), oar_radius*Y + oar_center(2), oar_radius*Z + oar_center(3));
set(h2, 'FaceAlpha', 0.8, 'FaceColor', 'magenta', 'EdgeColor', 'none');  % Magenta color for OAR
shading interp;

% Label OAR
text(oar_center(1), oar_center(2) + oar_radius, oar_center(3), 'OAR', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 8);

% Set plot limits based on extents of PTV and OAR
xlim([-150 150]);
ylim([-150 150]);
zlim([-150 150]);

% Plot Dose Box
box_center = dose_box(1:3);
box_side_lengths = dose_box(4:6);
O = box_center - box_side_lengths / 2;  % Origin of the box
plotcube(box_side_lengths, O, .1, [0.8 0.8 0]);  % Yellowish color for dose box

% Label Dose Box
text(box_center(1), box_center(2), box_center(3), 'Dose Box', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 8);

% Plot Isocenter
plot3(isocenter(1), isocenter(2), isocenter(3), '*r', 'MarkerSize', 10);  % Red star for isocenter

% Label Isocenter
text(isocenter(1), isocenter(2), isocenter(3), 'Isocenter', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 8);

% Plot Head
[Xhead, Yhead, Zhead] = ellipsoid(head_center(1), head_center(2), head_center(3), head_a, head_b, head_c, 20);
head = surf(Xhead, Yhead, Zhead);
set(head, 'FaceAlpha', 0.1, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none');  % Grey color for head

% Label Head
text(head_center(1), head_center(2) + head_b, head_center(3) + head_c, 'Head', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 8);

% Title for the plot
title('3-D Gamma Knife Radiosurgery Scene');

% Labels for axes
xlabel('x-axis (mm)');
ylabel('y-axis (mm)');
zlabel('z-axis (mm)');

hold off;

%% Compute Dose Absorption Function Table

fprintf('<strong>Part 2: Computing Dose Absorption Function Table</strong>\n')

% Determine the maximum head size
max_head_size = 2 * max([head_a, head_b, head_c]);
increment = 0.1; % mm

% Compute the dose absorption function table
dose_absorption_function_table = Compute_Depth_Dose(max_head_size, increment);

% Create the plot
figure();
plot(dose_absorption_function_table.Depth, dose_absorption_function_table.Dose, 'o-', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', [0, 0.4470, 0.7410]);
xlabel('Depth from Skin (mm)');
ylabel('Dose (%)');
title('Dose Absorption Function');
grid on;  % Enable grid for better readability

set(gca, 'FontSize', 10); % Set font size for axis labels
xlim([0, max_head_size]); % Set x-axis limits to the maximum head size
ylim([0, max(dose_absorption_function_table.Dose) * 1.1]); % Set y-axis limits slightly above the maximum dose

% Adding a line at depth 20 mm or 2 cm for reference
hold on;
line([20, 20], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
text(20, max(dose_absorption_function_table.Dose) * 0.9, 'Depth = 20 mm', 'Color', 'r', 'VerticalAlignment', 'bottom');

hold off;

%% Compute Radial Dose Function Table

fprintf('<strong>Part 3: Computing Radial Dose Absorption Function Table</strong>\n')

% Set parameters for the radial dose calculation
start_radial_dose = -30; % mm
end_radial_dose = 30; % mm
increment = 0.1; % mm

% Compute the radial dose function table
radial_dose_function_table = Compute_Radial_Dose(start_radial_dose, end_radial_dose, increment);

% Create the plot for the radial dose function
figure();
plot(radial_dose_function_table.Radial_Distance, radial_dose_function_table.Dose, 'o-', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', [0, 0.4470, 0.7410]);
xlabel('Radial Distance (mm)');
ylabel('Dose (%)');
title('Radial Dose Function');
grid on; % Enable grid for better readability

% Enhance plot appearance
set(gca, 'FontSize', 10); % Set font size for axis labels
xlim([start_radial_dose, end_radial_dose]); % Set x-axis limits to the range of radial distances
ylim([0, max(radial_dose_function_table.Dose) * 1.1]); % Set y-axis limits slightly above the maximum dose

hold off;

%% Compute Beam Directions

fprintf('<strong>Part 4: Computing Latitude, Longitude and unit direction vector for Beams</strong>\n')

beam_separation_angle = 30; %degrees
beam_diameter = 30; %mm
isocenter = [30, 0, 15]; %mm

[beam_nvectors] = Compute_Beam_Directions(beam_separation_angle);

% Latitudes and Longitudes calculated for for beam directions
longitudes = atan2(beam_nvectors(:, 2), beam_nvectors(:, 1));
latitudes = atan2(beam_nvectors(:, 3), sqrt(beam_nvectors(:, 1).^2 + beam_nvectors(:, 2).^2));

% Vectorized addition for origin points - make them arrays of the same length as beam_nvectors
u = ones(length(beam_nvectors), 1) * isocenter(1);
v = ones(length(beam_nvectors), 1) * isocenter(2);
w = ones(length(beam_nvectors), 1) * isocenter(3);

% Decompose beam vectors
[x, y, z] = deal(beam_nvectors(:, 1), beam_nvectors(:, 2), beam_nvectors(:, 3));

% First Quiver Plot
figure()
h7 = quiver3(u+x,v+y,w+z,-x,-y,-z);
set(h7,'AutoScale','on', 'AutoScaleFactor', 1.2)
xlabel('x')
ylabel('y')
zlabel('z')
title('Beam Directions Toward Isocenter')

% Second Quiver Plot
figure();
h7 = quiver3(u, v, w, x, y, z, 'b');
set(h7, 'AutoScale', 'on', 'AutoScaleFactor', 120);
xlabel('x');
ylabel('y');
zlabel('z');
title('Beam Directions');
hold on;

% Plot head
[Xhead, Yhead, Zhead] = ellipsoid(head_center(1), head_center(2), head_center(3), head_a, head_b, head_c, 20);
head = surf(Xhead, Yhead, Zhead);
set(head, 'FaceAlpha', 0.1);

%% Compute Skin Entry Points

fprintf('<strong>Part 5: Computing Skin Entry Points</strong>\n')

% Compute skin entry points
[skin_entry_points, depth_to_isocenter] = Compute_Skin_Entry_Points(beam_nvectors, isocenter, head_a, head_b, head_c);

figure();

% Plot skin entry points
plot3(skin_entry_points(:, 1), skin_entry_points(:, 2), skin_entry_points(:, 3), 'r*', 'MarkerSize', 10);
hold on;

% Plot head (ellipsoid)
[Xhead, Yhead, Zhead] = ellipsoid(head_center(1), head_center(2), head_center(3), head_a, head_b, head_c, 20);
head = surf(Xhead, Yhead, Zhead);
set(head, 'FaceAlpha', 0.1); % Set transparency for better visibility

% Plot beams
u = ones(size(beam_nvectors, 1), 1) * isocenter(1);
v = ones(size(beam_nvectors, 1), 1) * isocenter(2);
w = ones(size(beam_nvectors, 1), 1) * isocenter(3);
x = beam_nvectors(:, 1);
y = beam_nvectors(:, 2);
z = beam_nvectors(:, 3);

h8 = quiver3(u, v, w, x, y, z, 'b');
set(h8, 'AutoScale', 'on', 'AutoScaleFactor', 130);

% Enhance plot appearance
xlabel('x');
ylabel('y');
zlabel('z');
title('Skin Entry Points');
legend('Skin Entry Points', 'Location', 'best');
axis equal; % Set equal scaling for all axes
grid on; % Enable grid for better orientation

hold off;

%% Compute Beam Safety Flags

fprintf('<strong>Part 6: Computing Beam Safety Flags</strong>\n')

% Compute beam safety flags and determine unsafe beams
[beam_safety_flags, num_unsafe_beams] = Compute_Beam_Safety_Flags(oar_center, oar_radius, beam_diameter, isocenter, beam_nvectors);
indices_unsafe = find(beam_safety_flags);
unsafe_beams = beam_nvectors(indices_unsafe, :);

% Vectorized initialization for quiver plot origin
u = isocenter(1) * ones(size(unsafe_beams, 1), 1);
v = isocenter(2) * ones(size(unsafe_beams, 1), 1);
w = isocenter(3) * ones(size(unsafe_beams, 1), 1);

% Extract beam direction components
x = unsafe_beams(:, 1);
y = unsafe_beams(:, 2);
z = unsafe_beams(:, 3);

% Plot unsafe beams
figure()
h9 = quiver3(u,v,w,x,y,z,'b');
set(h9,'AutoScale','on', 'AutoScaleFactor', 120)
hold on

% Plot PTV
[X, Y, Z] = sphere(20);
h1 = surf(ptv_radius * X + ptv_center(1), ptv_radius * Y + ptv_center(2), ptv_radius * Z + ptv_center(3), 'FaceColor', 'green', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
text(ptv_center(1) + ptv_radius, ptv_center(2), ptv_center(3), 'PTV', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

% Plot OAR
h2 = surf(oar_radius * X + oar_center(1), oar_radius * Y + oar_center(2), oar_radius * Z + oar_center(3), 'FaceColor', 'blue', 'FaceAlpha', 0.8, 'EdgeColor', 'none');
text(oar_center(1), oar_center(2) + oar_radius, oar_center(3), 'OAR', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

% Enhance plot aesthetics
axis equal; % Equal scaling
xlim([-150 150]);
ylim([-150 150]);
zlim([-150 150]);
xlabel('x');
ylabel('y');
zlabel('z');
title('Beam Safety Flags');
hold off;

%% Compute Radial Distance

fprintf('<strong>Part 7: Computing Radial Distance</strong>\n')
fprintf('Testing the top helmet beam radial distance to lateral PTV edge:\n')

point_of_interest = [ptv_center(1)+ptv_radius,ptv_center(2),ptv_center(3)];
beam_index = 1;
radial_distance = Compute_Radial_Distance(point_of_interest, beam_nvectors, beam_index, isocenter)

%% Depth from Skin

fprintf('<strong>Part 8: Computing Depth from Skin</strong>\n')
fprintf('Testing the top helmet beam radial distance to lateral PTV edge:\n')

point_of_interest = [ptv_center(1)+ptv_radius,ptv_center(2),ptv_center(3)];
beam_index = 1;
skin_entry_point = skin_entry_points(beam_index,:)
[depth_from_skin] = Compute_Depth_from_Skin(point_of_interest, beam_nvectors, beam_index, isocenter, skin_entry_points)

%% Compute Point Dose from Beam

fprintf('<strong>Part 9: Computing Point Dose from Beam</strong>\n')
fprintf('Testing with the isocenter as the point of interest which is 60 mm in depth with an expected ~0.80 dose:\n')

[radial_distance] = Compute_Radial_Distance(isocenter, beam_nvectors, beam_index, isocenter);
[depth_from_skin] = Compute_Depth_from_Skin(isocenter, beam_nvectors, beam_index, isocenter, skin_entry_points);

[closest_depth,ind_closest_depth] = min(abs(dose_absorption_function_table.Depth-ones(length(dose_absorption_function_table.Depth),1).*depth_from_skin));
depth_dose = dose_absorption_function_table.Dose(ind_closest_depth)

[closest_rad,ind_closest_rad] = min(abs(radial_dose_function_table.Radial_Distance-ones(length(radial_dose_function_table.Radial_Distance),1).*radial_distance));
rad_dose = radial_dose_function_table.Dose(ind_closest_rad)

[point_dose_value] = Compute_Point_Dose_from_Beam(dose_absorption_function_table,depth_from_skin,radial_dose_function_table,radial_distance)

%% Compute Point Dose from All Beams
fprintf('<strong>Part 10A: Computing Point Dose from All Beams</strong>\n')
fprintf('Testing with the isocenter as the point of interest with all safe beams:\n')
point_of_interest_P13 = isocenter;

[full_isocenter_dose_value_from_safe_beams] = Compute_Dose_from_All_Beams(point_of_interest_P13,beam_nvectors,isocenter,skin_entry_points,dose_absorption_function_table,radial_dose_function_table,indices_unsafe)

fprintf('<strong>Part 10B: Computing Point Dose from All Beams to Check Dose Estimates</strong>\n')
fprintf('Testing with the isocenter as the point of interest with all beams:\n')
%The dose at the isocenter with ALL beams turned on
[isocenter_dose_value] = Compute_Dose_from_All_Beams(isocenter,beam_nvectors,isocenter,skin_entry_points,dose_absorption_function_table,radial_dose_function_table,[])

fprintf('Testing with the OAR point closest to PTV as the point of interest with all beams:\n')
%The dose at the point on the OAR closest to PTV with ALL beams turned on
offset_toward_isocenter = oar_radius*(oar_center-isocenter)/norm(oar_center-isocenter);
oar_closest_to_isocenter = oar_center+offset_toward_isocenter;
[OAR_max_dose_value] = Compute_Dose_from_All_Beams(oar_closest_to_isocenter,beam_nvectors,isocenter,skin_entry_points,dose_absorption_function_table,radial_dose_function_table,[])


%% Compute Surface Dose on PTV

fprintf('<strong>Part 11: Computing the Surface Dose to the PTV</strong>\n')

num_points = 2000;

[points_ptv,point_ptv_doses] = Compute_Surface_Dose(num_points,ptv_center,ptv_radius,beam_nvectors,isocenter,skin_entry_points,dose_absorption_function_table,radial_dose_function_table,indices_unsafe);

[max_ptv_surface_dose, ind_ptv_surface_dose] = max(point_ptv_doses);

%Create a colormap
[doseColors_ptv] = Create_Color_Map(points_ptv(:,1),point_ptv_doses);

figure()
scatter3(points_ptv(:,1), points_ptv(:,2), points_ptv(:,3), 14, doseColors_ptv, 'filled');
axis equal
colormap jet;
cb = colorbar;
cb.Title.String = "Dose";
clim([min(point_ptv_doses) max(point_ptv_doses)])

hold on
plot3(points_ptv(ind_ptv_surface_dose,1), points_ptv(ind_ptv_surface_dose,2), points_ptv(ind_ptv_surface_dose,3),'*r','linewidth', 10)
hold on
labelhighestdose = {['Highest Dose : ', num2str(max_ptv_surface_dose), 'units']};
text(points_ptv(ind_ptv_surface_dose,1), points_ptv(ind_ptv_surface_dose,2), points_ptv(ind_ptv_surface_dose,3), labelhighestdose,'VerticalAlignment','top','HorizontalAlignment','left')
title('PTV Surface Dose')
xlabel('x')
ylabel('y')
zlabel('z')

%% Compute Surface Dose on OAR

fprintf('<strong>Part 12: Computing the Surface Dose to the OAR</strong>\n')

num_points = 2000;

[points_oar,point_oar_doses] = Compute_Surface_Dose(num_points,oar_center,oar_radius,beam_nvectors,isocenter,skin_entry_points,dose_absorption_function_table,radial_dose_function_table,indices_unsafe);

[max_oar_surface_dose, ind_oar_surface_dose] = max(point_oar_doses);

%Create a colormap
[doseColors_oar] = Create_Color_Map(points_oar(:,1),point_oar_doses);
figure()
scatter3(points_oar(:,1), points_oar(:,2), points_oar(:,3), 14, doseColors_oar, 'filled');
axis equal
colormap jet;
cb = colorbar;
cb.Title.String = "Dose";
clim([min(point_oar_doses) max(point_oar_doses)])
hold on
plot3(points_oar(ind_oar_surface_dose,1), points_oar(ind_oar_surface_dose,2), points_oar(ind_oar_surface_dose,3),'*r','linewidth', 10)
hold on
xlabel('x')
ylabel('y')
zlabel('z')
labelhighestdose = {['Highest Dose : ', num2str(max_oar_surface_dose), 'units']};
text(points_oar(ind_oar_surface_dose,1), points_oar(ind_oar_surface_dose,2), points_oar(ind_oar_surface_dose,3), labelhighestdose,'VerticalAlignment','top','HorizontalAlignment','left')
title('OAR Surface Dose')

%% Compute Optimal Irradiation Time

% See document

%% Part 14A+15A - Dosimetry Analysis PTV and OAR
fprintf('<strong>Part 14A: Computing the Dose in the PTV</strong>\n')

voxel_size = DOSE_VOXEL_SIZE; %mm

%Build points of interest inside PTV box
[ptv_box] = Compute_Dose_Box(ptv_center, ptv_radius, ptv_center, ptv_radius);
[points_ptv,point_doses_ptv] = Compute_Dose(ptv_box, voxel_size, beam_nvectors,isocenter,skin_entry_points,dose_absorption_function_table,radial_dose_function_table,indices_unsafe);
[max_point_dose_ptv,index_max_point_dose_ptv] = max(point_doses_ptv);
[min_point_dose_ptv,~] = min(point_doses_ptv);

fprintf('<strong>Part 15A: Computing the Dose in the OAR</strong>\n')

%Build points of interest inside OAR box
[oar_box] = Compute_Dose_Box(oar_center, oar_radius, oar_center, oar_radius);
[points_oar,point_doses_oar] = Compute_Dose(oar_box, voxel_size, beam_nvectors,isocenter,skin_entry_points,dose_absorption_function_table,radial_dose_function_table,indices_unsafe);
[max_point_dose_oar,index_max_point_dose_oar] = max(point_doses_oar);

fprintf('<strong>Part 14A+15A: Plot the Dose Grid Boxes Around the PTV and the OAR</strong>\n')

%Create a colormaps
[doseColors] = Create_Color_Map([points_ptv(:,1);points_oar(:,1)],[point_doses_ptv;point_doses_oar]);

figure()
scatter3([points_ptv(:,1);points_oar(:,1)], [points_ptv(:,2);points_oar(:,2)], [points_ptv(:,3);points_oar(:,3)], 14, doseColors, 'filled');
axis equal
colormap jet
cb = colorbar;
cb.Title.String = "Dose";
clim([min([point_doses_ptv;point_doses_oar]) max([point_doses_ptv;point_doses_oar])])
xlabel('x')
ylabel('y')
zlabel('z')
title('Volume Dose in OAR and PTV Boxes')

%% Part 14B+15B - Dosimetry Analysis PTV and OAR Plots
fprintf('<strong>Part 14B+15B: Computing and Plotting the Dose Volume Histograms for PTV and OAR</strong>\n')

% Prescribed parameters
D100 = 20; % Units
DOARMAX = 10; % Maximum tolerable dose for OAR, adjust as necessary
dose_max = 29; % Units
inc = 0.01; % Units

% Creating Dose Volume Histogram for PTV
[relative_dose_ptv, ratio_of_total_structure_volume_ptv] = Create_Dose_Volume_Histogram(point_doses_ptv, D100, dose_max, inc);
figure();
plot(relative_dose_ptv, ratio_of_total_structure_volume_ptv, '-g');
xlabel('Relative Dose (%)');
ylabel('Ratio of Total Structure Volume (%)');
title('Dose Volume Histogram for PTV');
legend('PTV dose (relative to D100)');
xlim([0 150]);
grid on;

% Creating Dose Volume Histogram for OAR
[relative_dose_oar, ratio_of_total_structure_volume_oar] = Create_Dose_Volume_Histogram(point_doses_oar, DOARMAX, dose_max, inc);
figure();
plot(relative_dose_oar, ratio_of_total_structure_volume_oar, '-r');
xlabel('Relative Dose (%)');
ylabel('Ratio of Total Structure Volume (%)');
title('Dose Volume Histogram for OAR');
legend('OAR dose (relative to DOARMAX)');
xlim([0 150]);
grid on;

% Reporting over-dosage in OAR
fprintf('Number of OAR doses above DOARMAX:\n');
num_OAR_overdosed_points = sum(point_doses_oar > DOARMAX);
percentage_OAR_overdosed_points = num_OAR_overdosed_points / length(point_doses_oar);
fprintf('Count: %d, Percentage: %.2f%%\n', num_OAR_overdosed_points, percentage_OAR_overdosed_points * 100);