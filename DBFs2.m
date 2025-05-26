% Super-Resolution: MUSIC, ESPRIT vs Conventional Beamforming
clc, clear , close all
% Parameters
c							=	3e8;								% Speed of light (m/s)
f							=	10e9;								% Frequency (Hz)
lambda						=	c / f;								% Wavelength (m)
d							=	lambda / (1 + sind(45));			% Element spacing (half wavelength)
N							=	64;									% Number of array elements
angles						=	[20 20.8];							% True source angles (closely spaced < beamwidth)
SNR							=	15;									% High Signal-to-Noise Ratio (dB)
theta_range					=	-25:0.01:25;						% Angular range for analysis
snapshots					=	250;								% Increase number of snapshots
threshold_dB                =	-10;
min_distance                =   0.5;
		
%% Generate signals				
fs							=	1e6;								% Sampling frequency
t							=	(0:1/fs:(snapshots-1)/fs);			% Time vector
signal 						=	exp(1j * 2 * pi * f * t);			% Source 1
X							=	zeros(N, length(t));

% Add signals to array elements
for m                       =	1:length(angles)
	for n                   =	1:N
		X(n, :)             =	X(n, :) + signal * exp(-1j * 2 * pi * d * (n-1) * sind(-angles(m)) / lambda);
	end
end

% Add noise
X							=	awgn(X, SNR, 'measured');

% Compute the covariance matrix with diagonal loading
R							=	(X * X') / size(X, 2);				% Covariance matrix
R							=	R + eye(N) * 0.01;					% Diagonal loading for stability

%% Delay-and-Sum Beamforming
AF_DS						=	zeros(size(theta_range));
k							=	2 * pi / lambda;
for i						=	1:length(theta_range)
    theta					=	theta_range(i);
    steering_vector			=	exp(1j * k * d * (0:N-1)' * sind(theta));
    AF_DS(i)				=	abs(steering_vector' * mean(X, 2));	% Correct averaging
end
AF_DS						=	AF_DS / max(AF_DS);
det_DS                      =	merge_angles(theta_range, AF_DS, threshold_dB);
max_DS                      =	find_peaks(theta_range, AF_DS, threshold_dB, min_distance);

%% MVDR Algorithm
AF_MVDR						=	zeros(size(theta_range));
for i						=	1:length(theta_range)
    theta					=	theta_range(i);
    steering_vector			=	exp(1j * k * d * (0:N-1)' * sind(theta));
    AF_MVDR(i)				=	1 / real(steering_vector' / R * steering_vector);
end		
AF_MVDR						=	AF_MVDR / max(AF_MVDR);
det_MVDR					=	merge_angles(theta_range, AF_MVDR, threshold_dB);
max_MVDR					=	find_peaks(theta_range, AF_MVDR, threshold_dB, min_distance);

%% Beamscan Algorithm
AF_BS						=	zeros(size(theta_range));
for i						=	1:length(theta_range)
    theta					=	theta_range(i);
    steering_vector			=	exp(1j * k * d * (0:N-1)' * sind(theta));
    AF_BS(i)				=	abs(steering_vector' * R * steering_vector);
end
AF_BS						=	AF_BS / max(AF_BS);
det_BS                      =	merge_angles(theta_range, AF_BS, threshold_dB);
max_BS                      =	find_peaks(theta_range, AF_BS, threshold_dB, min_distance);

%% MUSIC Algorithm
[U, ~, ~]					=	svd(R);								% Eigen decomposition
Un							=	U(:, 3:end);						% Noise subspace
AF_MUSIC					=	zeros(size(theta_range));
for i						=	1:length(theta_range)
    theta					=	theta_range(i);
    steering_vector			=	exp(1j * k * d * (0:N-1)' * sind(theta));
    AF_MUSIC(i)				=	1 / abs(steering_vector' * Un * Un' * steering_vector);
end
AF_MUSIC					=	AF_MUSIC / max(AF_MUSIC);
det_MUSIC                   =	merge_angles(theta_range, AF_MUSIC, threshold_dB);
max_MUSIC                   =	find_peaks(theta_range, AF_MUSIC, threshold_dB, min_distance);
%% ESPRIT Algorithm
Es							=	U(:, 1:2);							% Signal subspace (for two sources)
phi							=	pinv(Es(1:end-1, :)) * Es(2:end, :);% Solve ESPRIT equation
eigenvalues					=	eig(phi);
det_ESPRIT                  =	asind(angle(eigenvalues) * lambda / (2 * pi * d));

%% Spatial Smoothing for MUSIC
R_smooth					=	zeros(8, 8);						% Initialize smoothed covariance matrix
L							=	N - 8 + 1;							% Number of valid subarrays
for l						=	1:L
    R_smooth				=	R_smooth + R(l:l+7, l:l+7);			% Average subarray covariance
end
R_smooth					=	R_smooth / L;						% Normalize by the number of subarrays

%% MUSIC with Spatial Smoothing
[U_smooth, ~, ~]			=	svd(R_smooth);						% Eigen decomposition
Un_smooth					=	U_smooth(:, 3:end);					% Noise subspace for smoothed matrix
AF_MUSIC_smooth				=	zeros(size(theta_range));
for i						=	1:length(theta_range)
    theta					=	theta_range(i);
    steering_vector			=	exp(1j * k * d * (0:7)' * sind(theta)); % Adjusted size for smoothing
    AF_MUSIC_smooth(i)		=	1 / abs(steering_vector' * Un_smooth * Un_smooth' * steering_vector);
end
AF_MUSIC_smooth				=	AF_MUSIC_smooth / max(AF_MUSIC_smooth);
det_MUSIC_smooth            =	merge_angles(theta_range, AF_MUSIC_smooth, threshold_dB);
max_MUSIC_smooth            =	find_peaks(theta_range, AF_MUSIC_smooth, threshold_dB, min_distance);
%% Plot results
figure;
plot(theta_range, 20*log10(AF_DS), 'LineWidth', 1.5); hold on;
plot(theta_range, 20*log10(AF_BS), 'LineWidth', 1.5);
plot(theta_range, 20*log10(AF_MUSIC), 'LineWidth', 1.5);
plot(theta_range, 20*log10(AF_MUSIC_smooth), '--', 'LineWidth', 1.5);
plot(theta_range, 20*log10(AF_MVDR), 'LineWidth', 1.5);
% xline(estimated_angles_ESPRIT, '--k', 'LineWidth', 1.5, 'Label', 'ESPRIT');
legend('Delay-and-Sum', 'Beamscan', 'MUSIC', 'MUSIC (Smoothed)', 'MVDR');%, 'ESPRIT');
xlabel('Angle (degrees)');
ylabel('Gain (dB)');
title('Super-Resolution vs Conventional Beamforming');
grid on;

disp(['True Angles are:', num2str(angles)]);
disp(['Angles detected by Delay and SUM are:', num2str(det_DS)])
disp(['Angles detected by MUSIC are:', num2str(det_MUSIC)])
disp(['Angles detected by Smooth MUSIC are:', num2str(det_MUSIC_smooth)])
disp(['Angles detected by MVDR are:', num2str(det_MVDR)])
disp(['Angles detected by Beamscan are:', num2str(det_BS)])
disp(['Angles detected by ESPRIT are:', num2str(det_ESPRIT')])

disp(['Maximum Peaks detected by Delay and SUM are:', num2str(max_DS.angles)])
disp(['Maximum Peaks detected by MUSIC are:', num2str(max_MUSIC.angles)])
disp(['Maximum Peaks detected by Smooth MUSIC are:', num2str(max_MUSIC_smooth.angles)])
disp(['Maximum Peaks detected by MVDR are:', num2str(max_MUSIC.angles)])
disp(['Maximum Peaks detected by Beamscan are:', num2str(max_BS.angles)])
