clear, clc, close all
set(0, 'DefaultFigureWindowStyle', 'docked')
phase_quantization									=	1;
%% %%%%%%%%%%%%%%%%%%%%%% SYSTEM PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fc													=	6e9;								% Central frequency (Hz)
c													=	3e8;								% Speed of light (m/s)
Azimuth_max											=	45;									% Maximum Azimuth Steering
Elevation_max										=	45;									% Maximum Elevation Steering
phase_shifter_ENOB									=	6;
SLL													=	-33;
nbar												=	4;
lambda												=	c/fc;								% Wavelength (m)
dx													=	lambda/(1 + sind(Azimuth_max));		% Element spacing within subarray (m)
dy													=	lambda/(1 + sind(Elevation_max));	% Element spacing within subarray (m)
k													=	2 * pi / lambda;					% Wave number
		
% Beam steering directions					
azimuth_steer										=	0;									% Azimuth angle to steer beam (degrees)
elevation_steer										=	0;									% Elevation angle to steer beam (degrees)

%% %%%%%%%%%%%%%%%%%%%%%% PSEUDO-CIRCULAR ARRAY CONFIGURATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aperture_distribution								=	[4, 4, 4, 4, 4, 4, 4, 4, 4, 4];		% Subarray Placements
total_rows											=	length(aperture_distribution);
max_subarrays										=	max(aperture_distribution);
				
subarray_size										=	[4, 4];							% [rows, columns] of elements per subarray
Tile_size											=	[4, 4];							% [rows, columns] of elements per Tile
elements_per_subarray								=	prod(subarray_size);
				
% Calculate total number of elements				
total_subarrays										=	sum(aperture_distribution);
total_elements										=	total_subarrays * elements_per_subarray;
N_element_AZ										=	max_subarrays*subarray_size(2);		% Number of elements in Azimuth Coordinate
N_element_EL										=	total_rows*subarray_size(1);		% Number of elements in Elevation Coordinate
fprintf('Array Configuration:\n');
fprintf('  Total subarrays: %d\n', total_subarrays);
fprintf('  Elements per subarray: %d\n', elements_per_subarray);
fprintf('  Total elements: %d\n', total_elements);

%% %%%%%%%%%%%%%%%%%%%%%% CREATE ARRAY GEOMETRY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
element_positions									=	zeros(total_elements,2);
subarray_centers									=	zeros(total_subarrays,2);
element_count										=	0;
iii                                 				=   0;
jjj                                 				=   0;
max_x = 0;		
max_y = 0;		
for row												=	1:total_rows
    num_subarrays_in_row							=	aperture_distribution(row);
    
    % column positions for this row (centered)
    col_offset										=	-(num_subarrays_in_row - 1) * subarray_size(2) * dx / 2;
					
    for col											=	1:num_subarrays_in_row
        iii                         				=   iii + 1;
				
		% Subarray center position				
        subarray_x									=	col_offset + (col-1) * subarray_size(2) * dx;
        subarray_y									=	(row-1) * subarray_size(1) * dy;
        subarray_centers(iii,:)     				=   [subarray_x, subarray_y];
        
		% Create elements within this subarray
        for i										=	1:subarray_size(1)
            for j									=	1:subarray_size(2)
                element_x							=	subarray_x + (j - (subarray_size(2)+1)/2) * dx;
                element_y							=	subarray_y + (i - (subarray_size(1)+1)/2) * dy;
                jjj                 				=   jjj + 1;
                element_positions(jjj,:)			=	[element_x, element_y];
                element_count						=	element_count + 1;
																	
            end		
        end		
    end		
end		
			
% Center the entire array		
element_positions(:,1)								=	element_positions(:,1) - mean(element_positions(:,1));
element_positions(:,2)								=	element_positions(:,2) - mean(element_positions(:,2));
subarray_centers(:,1)								=   subarray_centers(:,1) - mean(subarray_centers(:,1));
subarray_centers(:,2)								=   subarray_centers(:,2) - mean(subarray_centers(:,2));
		
Tile_centers										=	subarray_centers;


%% %%%%%%%%%%%%%%%%%%%%%% PLOT ARRAY GEOMETRY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scatter(element_positions(:,1)/lambda, element_positions(:,2)/lambda, 20, 'filled');
axis equal; grid on; hold on;
scatter(subarray_centers(:,1)/lambda, subarray_centers(:,2)/lambda, 100, 'ro', 'LineWidth', 2);
scatter(Tile_centers(:,1)/lambda, Tile_centers(:,2)/lambda, 50, 'ko', 'LineWidth', 1);
xlabel('X Position (wavelengths)'); ylabel('Y Position (wavelengths)'); legend('Elements', 'Subarray Centers', 'Tile Centers');
title(sprintf('Pseudo-Circular Array Geometry\n%d Subarrays, %d Total Elements', total_subarrays, total_elements));

%% %%%%%%%%%%%%%%%%%%%%%% BEAMFORMING CALCULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
az_rad												=	deg2rad(azimuth_steer);
el_rad												=	deg2rad(elevation_steer);
				
u_steer												=	sin(az_rad);%cos(el_rad) * cos(az_rad);				
v_steer												= 	sin(el_rad);%cos(el_rad) * sin(az_rad);				
		
% Calculate steering vector (phase shifts)		
phase_shifts										=	k * (element_positions(:,1) * u_steer + element_positions(:,2) * v_steer);
if phase_quantization								==	1
	phase_steps										=	2*pi/(2^phase_shifter_ENOB);
	phase_shifts									=	round(phase_shifts/phase_steps)*phase_steps;
end		
steering_vector										=	exp(1j * phase_shifts);
weights												=	steering_vector / sqrt(total_elements);

%% %%%%%%%%%%%%%%%%%%%%%% Taper Design %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_elements_x										=	N_element_AZ;%/Tile_size(2);
max_elements_y										=	N_element_EL;%/Tile_size(1);
		
% Generate 1D Taylor and Bayliss windows		
taylor_x_tmp										=	taylorwin(N_element_AZ/Tile_size(2), nbar, SLL);
taylor_y_tmp										=	taylorwin(N_element_EL/Tile_size(1), nbar, SLL);
bayliss_x_tmp										=	baylisswin(N_element_AZ/Tile_size(2), nbar, SLL);
bayliss_y_tmp										=	baylisswin(N_element_EL/Tile_size(1), nbar, SLL);
for i 												=	1:length(taylor_x_tmp)
	taylor_x((i-1)*Tile_size(2)+1:i*Tile_size(2))	=	repmat(taylor_x_tmp(i),Tile_size(2),1);
	bayliss_x((i-1)*Tile_size(2)+1:i*Tile_size(2))	=	repmat(bayliss_x_tmp(i),Tile_size(2),1);
end
for i 												=	1:length(taylor_y_tmp)
	taylor_y((i-1)*Tile_size(1)+1:i*Tile_size(1))	=	repmat(taylor_y_tmp(i),Tile_size(1),1);
	bayliss_y((i-1)*Tile_size(1)+1:i*Tile_size(1))	=	repmat(bayliss_y_tmp(i),Tile_size(1),1);
end
% Map element positions to indices
x_positions											=	element_positions(:,1);
y_positions											=	element_positions(:,2);
x_min												=	min(x_positions); 
x_max												=	max(x_positions);
y_min												=	min(y_positions); 
y_max												=	max(y_positions);
		
if x_max											==	x_min
    x_indices										=	ones(size(x_positions)) * (max_elements_x/2);
else										
    x_indices										=	round((x_positions - x_min) / (x_max - x_min) * (max_elements_x-1)) + 1;
end		
		
if y_max											==	y_min
    y_indices										=	ones(size(y_positions)) * (max_elements_y/2);
else										
    y_indices										=	round((y_positions - y_min) / (y_max - y_min) * (max_elements_y-1)) + 1;
end						   		
		
% Initialize weights for three cases		
weights_tt											=	zeros(total_elements, 1);
weights_tb											=	zeros(total_elements, 1);
weights_bt											=	zeros(total_elements, 1);
		
% Assign weights		
for i												=	1:total_elements
    weights_tt(i)									=	taylor_x(x_indices(i)) * taylor_y(y_indices(i));
    weights_tb(i)									=	taylor_x(x_indices(i)) * bayliss_y(y_indices(i));
    weights_bt(i)									=	bayliss_x(x_indices(i)) * taylor_y(y_indices(i));
end		
		
% Combine with steering vector and normalize		
weights_tt											=	steering_vector .* weights_tt;
weights_tb											=	steering_vector .* weights_tb;
weights_bt											=	steering_vector .* weights_bt;
												
weights_tt											=	weights_tt / sqrt(sum(abs(weights_tt).^2));
weights_tb											=	weights_tb / sqrt(sum(abs(weights_tb).^2));
weights_bt											=	weights_bt / sqrt(sum(abs(weights_bt).^2));

%% %%%%%%%%%%%%%%%%%%%%%% PATTERN CALCULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
az_angles											=	-60:0.1:60;                 
el_angles											=	-60:0.1:60;                 
pattern_tt											=	zeros(length(az_angles), length(el_angles));
pattern_tb											=	zeros(length(az_angles), length(el_angles));
pattern_bt											=	zeros(length(az_angles), length(el_angles));
				
for az_idx											=	1:length(az_angles)
    for el_idx										=	1:length(el_angles)
        az_rad_current								=	deg2rad(az_angles(az_idx));
        el_rad_current								=	deg2rad(el_angles(el_idx));
        u_current									=	sin(az_rad_current);
        v_current									=	sin(el_rad_current);
																
			
        array_response								=	exp(1j * k * (element_positions(:,1) * u_current + element_positions(:,2) * v_current));
        pattern_tt(az_idx, el_idx)					=	abs(weights_tt' * array_response);
        pattern_tb(az_idx, el_idx)					= 	abs(weights_tb' * array_response);
        pattern_bt(az_idx, el_idx)					=	abs(weights_bt' * array_response);
    end
end

%% %%%%%%%%%%%%%%%%%%%%%% PLOT PATTERNS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Taylor × Taylor (sum pattern)
[~, el_idx]											=	min(abs(el_angles - elevation_steer));
az_pattern_tt										=	pattern_tt(:, el_idx);
figure; plot(az_angles, 20*log10(az_pattern_tt)); grid on;
xlabel('Azimuth Angle (degrees)'); ylabel('Gain (dB)');
title(sprintf('Azimuth Pattern at Elevation = %.1f°', elevation_steer)); xlim([-90, 90]); ylim([-50,50]);

[~, az_idx]											=	min(abs(az_angles - azimuth_steer));
el_pattern_tt										=	pattern_tt(az_idx, :);
figure; plot(el_angles, 20*log10(el_pattern_tt)); grid on;
xlabel('Elevation Angle (degrees)'); ylabel('Gain (dB)');
title(sprintf('Elevation Pattern at Azimuth = %.1f°', azimuth_steer)); xlim([-90, 90]); ylim([-50,50]);

% Taylor × Bayliss (difference in elevation)
el_pattern_tb										=	pattern_tb(az_idx, :);											  
figure; plot(az_angles, 20*log10(el_pattern_tb)); grid on;
xlabel('Azimuth Angle (degrees)'); ylabel('Gain (dB)');
title(sprintf('Elevation Difference Pattern at Azimuth = %.1f°', azimuth_steer)); xlim([-90, 90]); ylim([-50,50]);																  

% Bayliss × Taylor (difference in azimuth)
az_pattern_bt										=	pattern_bt(:, el_idx);																				   
figure; plot(el_angles, 20*log10(az_pattern_bt)); grid on;
xlabel('Elevation Angle (degrees)'); ylabel('Gain (dB)');
title(sprintf('Azimuth Difference Pattern at Elevation = %.1f°', elevation_steer)); xlim([-90, 90]); ylim([-50,50]);

% 3D Pattern
figure; [AZ, EL]									=	meshgrid(az_angles, el_angles);
mesh(AZ, EL, 20*log10(pattern_tt'));
xlabel('Azimuth (degrees)'); ylabel('Elevation (degrees)'); zlabel('Gain (dB)');
title('3D Radiation SUM Pattern'); colorbar; view(45, 30); axis tight;

figure; [AZ, EL]									=	meshgrid(az_angles, el_angles);
mesh(AZ, EL, 20*log10(pattern_tb'));
xlabel('Azimuth (degrees)'); ylabel('Elevation (degrees)'); zlabel('Gain (dB)');
title('3D Radiation Elevation Difference Pattern'); colorbar; view(45, 30); axis tight;

figure; [AZ, EL]									=	meshgrid(az_angles, el_angles);
mesh(AZ, EL, 20*log10(pattern_bt'));
xlabel('Azimuth (degrees)'); ylabel('Elevation (degrees)'); zlabel('Gain (dB)');
title('3D Radiation Azimuth Difference Pattern'); colorbar; view(45, 30); axis tight;

%% %%%%%%%%%%%%%%%%%%%%%% BEAMPOINTING CALCULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[peak_az_idx,peak_el_idx]							=	find(pattern_tt == max(pattern_tt(:)));
actual_az_steer										=	az_angles(peak_az_idx);
actual_el_steer										=	el_angles(peak_el_idx);
fprintf('Requested steering: Az=%.1f°, El=%.1f°\n', azimuth_steer, elevation_steer);
fprintf('Actual steering: Az=%.1f°, El=%.1f°\n', actual_az_steer, actual_el_steer);

%% %%%%%%%%%%%%%%%%%%%%%% BEAMWIDTH CALCULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate 3dB beamwidth for azimuth pattern
az_pattern_db										=	20*log10(az_pattern_tt);
peak_gain											=	max(az_pattern_db);
half_power											=	peak_gain - 3;
				
% Find -3dB points				
above_3db											=	az_pattern_db >= half_power;
indices												=	find(diff(above_3db) ~= 0);
az_bw												=	az_angles(indices(2)) - az_angles(indices(1));
fprintf('Azimuth Beamwidth: %.2f degrees\n', az_bw);	

% Calculate 3dB beamwidth for elevation pattern
el_pattern_db										=	20*log10(el_pattern_tt);
peak_gain											=	max(el_pattern_db);
half_power											=	peak_gain - 3;
				
% Find -3dB points				
above_3db											=	el_pattern_db >= half_power;
indices												=	find(diff(above_3db) ~= 0);
el_bw												=	el_angles(indices(2)) - el_angles(indices(1));
fprintf('Elevation Beamwidth: %.2f degrees\n', el_bw);
fprintf('Peak Gain: %.2f dB\n', peak_gain);

function w											=	baylisswin(N, nbar, sll) 
	sll 											=	sll - 20;
    n                                       		=   (1:N)' - (N+1)/2;
	w_t                                     		=   taylorwin(N,nbar,sll);
    w                                       		=   w_t.*n;
    w												=	w/sqrt(sum(abs(w).^2));
    % A												=	acosh(10^(-sll/20))/pi;
	% w												=	zeros(N,1);
	% sigma											=	1.2;
	% xi											=	sinh(A)/(N/2);
	% u_k											=	zeros(nbar,1);
	% for k											=	1:nbar
	% 	u_k(k)										=	sigma*xi*(1 + (k-0.5)/(nbar+1));
	% end		
    % for k											=	1:nbar
    %     w											=	w + (-1)^(k+1)*sin(pi + u_k(k)*(2*n+1)/N);
    % end		
	
    % w												=   w_t.*w;
	% w												=	w.*sign(sum(abs(w).^2));
end