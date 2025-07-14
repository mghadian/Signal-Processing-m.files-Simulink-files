% INS implementation = IMU + processor + software algorithms (+ sometimes GPS)
% Direct Implementation of GPS/INS navigation using ESKF (Error-state Kalman Filter)
clear ; close all; clc;

%% ====================== LOAD DATA ======================
% Replace with your actual data loading
% Assume we have the following variables:
% imu_data: [time, acc_x, acc_y, acc_z, gyro_x, gyro_y, gyro_z]
% gps_data: [time, lat, lon, alt, v_n, v_e, v_d]

fs_imu										=	100;				% IMU sample rate (Hz)
fs_gps										= 	1;					% GPS sample rate (Hz)
duration									=	300;				% seconds
t_imu										=	(0:1/fs_imu:duration)';
t_gps										=	(0:1/fs_gps:duration)';

%% Vehicle dynamic scenario
N											=	length(t_imu);

% Circular trajectory (constant velocity)
omega										=	0.5;				% Angular velocity (rad/s)
radius										=	10;					% Circle radius (m)

% Position, velocity, acceleration
pos_gt										=	radius * [cos(omega * t_imu), sin(omega * t_imu), zeros(N,1)];
vel_gt										=	radius * omega * [-sin(omega * t_imu), cos(omega * t_imu), zeros(N,1)];
acc_gt										=	radius * omega^2 * [-cos(omega * t_imu), -sin(omega * t_imu), zeros(N,1)];

%% Generate synthetic IMU data (acc in m/s², gyro in rad/s)
acc_gt										=	acc_gt + [0, 0, 9.81];	% Include gravity

% Simulate IMU noise and biases
acc_bias_true								=	[0.1, -0.05, 0.2];		% True constant bias (m/s²)
gyro_bias_true								=	[0.01, -0.01, 0.005];	% True gyro bias (rad/s)

% Generate noisy IMU measurements
acc_meas									=	acc_gt + acc_bias_true + 0.05 * randn(N,3);
gyro_meas									=	[zeros(N,2), omega * ones(N,1)] + gyro_bias_true + 0.01 * randn(N,3);

% Combine into IMU data matrix
imu_data									=	[t_imu, acc_meas, gyro_meas];

%% Generate synthetic GPS data (lat/lon in deg, alt in m, vel in m/s)
M											=	length(t_gps);

% Sample ground truth at GPS timestamps
[~, idx]									=	ismember(t_gps, t_imu);
pos_gps_gt									=	pos_gt(idx, :);
vel_gps_gt									=	vel_gt(idx, :);

% Add GPS noise
pos_noise									=	1.0 * randn(M,3);		% Position noise (1m std)
vel_noise									=	0.1 * randn(M,3);		% Velocity noise (0.1 m/s std)

% Noisy GPS measurements
gps_pos_meas								=	pos_gps_gt + pos_noise;
gps_vel_meas								=	vel_gps_gt + vel_noise;

% Combine into GPS data matrix
gps_data									=	[t_gps, gps_pos_meas, gps_vel_meas];

%% ====================== INITIALIZATION ======================
% Initial state (position in NED, velocity in NED, attitude quaternion)
p											=	[0; 0; 0];          % NED position (m)
v											=	[0; 0; 0];          % NED velocity (m/s)
q											=	[1; 0; 0; 0];       % Quaternion (body to NED)

% Initial bias estimates (start with zero)
acc_bias									=	[0; 0; 0];			% Accelerometer bias (m/s²)
gyro_bias									=	[0; 0; 0];			% Gyro bias (rad/s)

% Earth parameters
omega_ie									=	7.292115e-5;		% Earth rotation rate (rad/s)(omega_ie*24*60*60/2pi = 1)
R0											=	6378137;			% Earth radius (m)
e											=	0.08181919;         % Earth eccentricity
g											=	9.81;               % Gravity (m/s²)
% Reference position for NED conversion
ref_lla										=	[0; 0; 0];			% Reference at origin for simplicity	
								   
% Kalman filter initialization
dt_imu										=	1/fs_imu;			% IMU time step

% State vector
x_est										=	zeros(15,1);		% Error state estimate
P											=	eye(15);            % Error covariance matrix

% Process noise (Q) and Measurement noise (R)
Q											=	diag([0.001^2, 0.001^2, 0.001^2, ...				% Position
												0.01^2, 0.01^2, 0.01^2, ...							% Velocity
												(0.5*pi/180)^2, (0.5*pi/180)^2, (0.5*pi/180)^2, ... % Attitude
												(1e-5)^2, (1e-5)^2, (1e-5)^2, ...					% Acc bias
												(1e-6)^2, (1e-6)^2, (1e-6)^2]);						% Gyro bias

R											=	diag([2^2, 2^2, 2^2, ...							% GPS position (m)
												0.2^2, 0.2^2, 0.2^2]);								% GPS velocity (m/s)

% Transformation matrices
C_bn										=	quat2dcm(q');                       				% Body to NED DCM (Direction Cosine Matrix)

% Preallocate results for plotting
results.time								=	t_imu;
results.p_ins								=	zeros(length(t_imu),3);
results.v_ins								=	zeros(length(t_imu),3);
results.att_ins								=	zeros(length(t_imu),3);
results.p_gps								=	zeros(length(t_imu),3);
results.v_gps								=	zeros(length(t_imu),3);
results.p_gt								=	pos_gt;

%% ====================== MAIN LOOP ======================
gps_idx										=	1;
for k										=	1:length(t_imu)
    % Current IMU measurement
    acc										=	imu_data(k,2:4)' - acc_bias;
    gyro									=	imu_data(k,5:7)' - gyro_bias;
    
    % ========== INS MECHANIZATION ==========
    % Attitude update
    omega									=	gyro;
    q										=	quat_update(q, omega, dt_imu);
	q										=	q / norm(q);
    C_bn									=	quat2dcm(q');
    
    % Velocity update (in NED frame)
    acc_n									=	C_bn * acc;
    acc_n									=	acc_n - [0; 0; g];  % Remove gravity
    v										=	v + acc_n * dt_imu;
    
    % Position update
    p										=	p + v * dt_imu;
    
    % Store results
    results.p_ins(k,:)						=	p';
    results.v_ins(k,:)						=	v';
    euler									=	quat2eul_manual(q);
    results.att_ins(k,:)					=	euler' * 180/pi;
	
	% ========== KALMAN FILTER PREDICTION ==========
    % State transition matrix
    F										=	build_F_matrix(C_bn, acc, dt_imu);
    
    % Predict error state and covariance
    x_est									=	F * x_est;
    P										=	F * P * F' + Q * dt_imu;												
    % ========== KALMAN FILTER UPDATE ==========
    % Only update when GPS measurement is available
    if gps_idx								<=	length(t_gps) && abs(t_imu(k) - t_gps(gps_idx)) < dt_imu/2
		% GPS measurement
		gps_p								=	gps_data(gps_idx,2:4).';
		gps_v								=	gps_data(gps_idx,5:7).';
		
		% Measurement innovation
		z									= [gps_p - p; 
												gps_v - v];
		
		% Proper measurement matrix
		H									=	[eye(3), zeros(3,12);           % Position observation
												zeros(3), eye(3), zeros(3,9)];	% Velocity observation
        
		% Kalman update
		S									=	H * P * H' + R;
        K									=	P * H' / S; 
		x_est								=	x_est + K * (z - H * x_est);
		P									=	(eye(15) - K * H) * P;
		
		% Apply corrections
		p									=	p + x_est(1:3);
		v									=	v + x_est(4:6);
        delta_theta							=	x_est(7:9);
        if norm(delta_theta)				>	0
            delta_q							=	[cos(norm(delta_theta)/2); 
												sin(norm(delta_theta)/2) * delta_theta/norm(delta_theta)];
            q								=	quatmultiply_manual(q, delta_q);
            q								=	q / norm(q);
            C_bn							=	quat2dcm(q');
        end
		
		% Bias correction
		acc_bias							=	acc_bias + x_est(10:12);
		gyro_bias							=	gyro_bias + x_est(13:15);
		
        % Reset error state after correction
        x_est                               =	zeros(15,1);										  
		% Store GPS data
		results.p_gps(k,:)					=	gps_p';
		results.v_gps(k,:)					=	gps_v';
		
		gps_idx								=	gps_idx + 1;
    else
        % Store GPS data (interpolated or held)
        if gps_idx							>	1
            results.p_gps(k,:)				=	results.p_gps(k-1,:);
            results.v_gps(k,:)				=	results.v_gps(k-1,:);
        end
    end

end

%% ====================== PLOTTING ======================
figure('Name','Position Estimation');
subplot(311); 
plot(results.time, results.p_gt(:,1), 'k--', 'LineWidth', 1.5); hold on;
plot(results.time, results.p_ins(:,1), 'b', 'LineWidth', 1);
plot(results.time, results.p_gps(:,1), 'r.', 'MarkerSize', 8);
title('North Position'); xlabel('Time (s)'); ylabel('m'); 
legend('Ground Truth', 'INS', 'GPS', 'Location', 'best');
grid on;

subplot(312); 
plot(results.time, results.p_gt(:,2), 'k--', 'LineWidth', 1.5); hold on;
plot(results.time, results.p_ins(:,2), 'b', 'LineWidth', 1);
plot(results.time, results.p_gps(:,2), 'r.', 'MarkerSize', 8);
title('East Position'); xlabel('Time (s)'); ylabel('m');
grid on;

subplot(313); 
plot(results.time, results.p_gt(:,3), 'k--', 'LineWidth', 1.5); hold on;
plot(results.time, results.p_ins(:,3), 'b', 'LineWidth', 1);
plot(results.time, results.p_gps(:,3), 'r.', 'MarkerSize', 8);
title('Down Position'); xlabel('Time (s)'); ylabel('m');
grid on;

figure('Name','Velocity Estimation');
subplot(311); 
plot(results.time, vel_gt(:,1), 'k--', 'LineWidth', 1.5); hold on;
plot(results.time, results.v_ins(:,1), 'b', 'LineWidth', 1);
plot(results.time, results.v_gps(:,1), 'r.', 'MarkerSize', 8);
title('North Velocity'); xlabel('Time (s)'); ylabel('m/s'); 
legend('Ground Truth', 'INS', 'GPS', 'Location', 'best');
grid on;

subplot(312); 
plot(results.time, vel_gt(:,2), 'k--', 'LineWidth', 1.5); hold on;
plot(results.time, results.v_ins(:,2), 'b', 'LineWidth', 1);
plot(results.time, results.v_gps(:,2), 'r.', 'MarkerSize', 8);
title('East Velocity'); xlabel('Time (s)'); ylabel('m/s');
grid on;

subplot(313); 
plot(results.time, vel_gt(:,3), 'k--', 'LineWidth', 1.5); hold on;
plot(results.time, results.v_ins(:,3), 'b', 'LineWidth', 1);
plot(results.time, results.v_gps(:,3), 'r.', 'MarkerSize', 8);
title('Down Velocity'); xlabel('Time (s)'); ylabel('m/s');
grid on;

figure('Name','Attitude Estimation');
subplot(311); plot(results.time, results.att_ins(:,1));
title('Yaw Angle'); xlabel('Time (s)'); ylabel('deg'); grid on;
subplot(312); plot(results.time, results.att_ins(:,2));
title('Pitch Angle'); xlabel('Time (s)'); ylabel('deg'); grid on;
subplot(313); plot(results.time, results.att_ins(:,3));
title('Roll Angle'); xlabel('Time (s)'); ylabel('deg'); grid on;

%% ====================== ERROR ANALYSIS ======================
% Calculate RMS errors against ground truth
pos_error_gt								= results.p_ins - results.p_gt;
vel_error_gt								= results.v_ins - vel_gt;

% Calculate RMS errors against GPS (only when GPS is available)
valid_idx									=	find(results.p_gps(:,1) ~= 0);
pos_error_gps								=	results.p_ins(valid_idx,:) - results.p_gps(valid_idx,:);
vel_error_gps								=	results.v_ins(valid_idx,:) - results.v_gps(valid_idx,:);

rms_pos_gt									=	sqrt(mean(pos_error_gt.^2));
rms_vel_gt									=	sqrt(mean(vel_error_gt.^2));
rms_pos_gps									=	sqrt(mean(pos_error_gps.^2));
rms_vel_gps									=	sqrt(mean(vel_error_gps.^2));

fprintf('=== RMS Errors vs Ground Truth ===\n');
fprintf('Position RMS Errors:\n North: %.3f m, East: %.3f m, Down: %.3f m\n', rms_pos_gt);
fprintf('Velocity RMS Errors:\n North: %.3f m/s, East: %.3f m/s, Down: %.3f m/s\n', rms_vel_gt);

fprintf('\n=== RMS Errors vs GPS ===\n');
fprintf('Position RMS Errors:\n North: %.3f m, East: %.3f m, Down: %.3f m\n', rms_pos_gps);
fprintf('Velocity RMS Errors:\n North: %.3f m/s, East: %.3f m/s, Down: %.3f m/s\n', rms_vel_gps);

fprintf('\n=== Final Bias Estimates ===\n');
fprintf('Accelerometer bias: [%.4f, %.4f, %.4f] m/s²\n', acc_bias);
fprintf('Gyroscope bias: [%.6f, %.6f, %.6f] rad/s\n', gyro_bias);
fprintf('True acc bias: [%.4f, %.4f, %.4f] m/s²\n', acc_bias_true);
fprintf('True gyro bias: [%.6f, %.6f, %.6f] rad/s\n', gyro_bias_true);

%% ====================== HELPER FUNCTIONS ======================
function q_new								=	quat_update(q, omega, dt)
    % Quaternion update using angular rate
    omega_norm								=	norm(omega);
    if omega_norm							>	eps
        axis								=	omega/omega_norm;
        angle								=	omega_norm*dt;
        dq									=	[cos(angle/2); axis*sin(angle/2)];
        q_new								=	quatmultiply_manual(q, dq);
    else
        q_new								=	q;
    end
    q_new									=	q_new / norm(q_new);
end

function F									=	build_F_matrix(C_bn, acc, dt)
    % Build state transition matrix for ESKF
    F										=	eye(15);
    
    % Position-velocity coupling
    F(1:3,4:6)								=	eye(3)*dt;
    
    % Velocity-attitude coupling
    F(4:6,7:9)								=	-skew_symmetric(C_bn*acc)*dt;
    
    % Velocity-accelerometer bias coupling
    F(4:6,10:12)							=	-C_bn*dt;
    
    % Attitude-gyroscope bias coupling
    F(7:9,13:15)							=	-C_bn*dt;

end

function S									=	skew_symmetric(v)
    % Skew-symmetric matrix
    S										=	[0, -v(3), v(2);
												v(3), 0, -v(1);
												-v(2), v(1), 0];
end

function euler								=	quat2eul_manual(q)
    % Manual quaternion to Euler angles conversion (ZYX order)
    q0										=	q(1); 
	q1										=	q(2); 
	q2										=	q(3); 
	q3										=	q(4);
    
    % Roll (x-axis rotation)
    sinr_cosp								=	2*(q0*q1 + q2*q3);
    cosr_cosp								=	1 - 2*(q1^2 + q2^2);
    roll									=	atan2(sinr_cosp, cosr_cosp);
    
    % Pitch (y-axis rotation)
    sinp									=	2*(q0*q2 - q3*q1);
    if abs(sinp)							>=	1
        pitch								=	sign(sinp)*pi/2;
    else
        pitch								=	asin(sinp);
    end
    
    % Yaw (z-axis rotation)
    siny_cosp								=	2*(q0*q3 + q1*q2);
    cosy_cosp								=	1 - 2*(q2^2 + q3^2);
    yaw										=	atan2(siny_cosp, cosy_cosp);
    
    euler									=	[yaw; pitch; roll]; % ZYX order
end

function q									=	quatmultiply_manual(q1, q2)
    % Manual quaternion multiplication
    q										=	[q1(1)*q2(1) - q1(2)*q2(2) - q1(3)*q2(3) - q1(4)*q2(4);
												q1(1)*q2(2) + q1(2)*q2(1) + q1(3)*q2(4) - q1(4)*q2(3);
												q1(1)*q2(3) - q1(2)*q2(4) + q1(3)*q2(1) + q1(4)*q2(2);
												q1(1)*q2(4) + q1(2)*q2(3) - q1(3)*q2(2) + q1(4)*q2(1)];
end
