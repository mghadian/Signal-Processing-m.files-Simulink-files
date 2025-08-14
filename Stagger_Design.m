% Design filter coefficients for staggered PRF MTI canceler using least squares optimization with adjustable weight W.
clc, close all

%% Parameters
PRF_avg									=	4e3/1.075;											% Maximum PRF in Hz
f_min									=	680;												% Min Doppler frequency in Hz (start of passband)
f_max									=	34000;												% Max Doppler frequency in Hz (end of passband)
f_c										=	f_min * 0.8;										% Stopband cutoff (80% of f_min, adjustable)
target_avg_dB							=	-45;												% Target average stopband magnitude in dB			

%% Normalize frequencies						
f_stop_norm								=	linspace(0, f_c / PRF_avg, 200);					% Discretized stopband
f_pass_norm								=	linspace(f_min / PRF_avg, f_max / PRF_avg, 1000);	% Discretized passband

%% Example stagger sequence lengths and ratios
stagger_examples						=	{
											struct('M', 2, 'ratios', [8, 9]), ... 				% Sequence length 2 (3 taps)
											struct('M', 3, 'ratios', [8, 9, 10]), ... 			% Sequence length 3 (4 taps) (or 24,25,27)
											struct('M', 4, 'ratios', [25, 30, 27, 31]), ... 	% Sequence length 4 (5 taps)
											struct('M', 5, 'ratios', [51, 62, 52, 60, 55]) ... 	% Sequence length 5 (6 taps)
											};
		
for ex									=	1:length(stagger_examples)
    M									=	stagger_examples{ex}.M;
    n									=	stagger_examples{ex}.ratios;
    N									=	M + 1;												% Number of taps = M + 1
		
    % Compute PRI values		
    sum_n								=	sum(n);
    average_n							=	sum_n / M;
    T_avg								=	1 / PRF_avg;
    T_basic								=	T_avg / average_n;
    T									=	n * T_basic;
    tau									=	cumsum([0, T]);										% Cumulative time delays
    tau_norm							=	tau / T_avg;										% Normalized

    % First blind frequency
    g									=	n(1);
    for k								=	2:M
        g								=	gcd(g, n(k));
    end		
    first_blind_factor					=	sum_n / (M * g);
    first_blind							=	first_blind_factor * PRF_avg;

    % Matrices for stopband and passband (complex)
    Astop								= 	exp(-1j * 2 * pi * f_stop_norm' * tau_norm);
    Apass								= 	exp(-1j * 2 * pi * f_pass_norm' * tau_norm);

    % Constraints for second-order null (K=2): sum x = 0, sum tau x = 0
    Aeq									=	[ones(1, N); tau_norm];
    beq									=	[0; 0];
    lb									=	-inf * ones(N, 1);
    ub									=	inf * ones(N, 1);

    % Bisection search for W to achieve target average dB
    W_low								=	1;
    W_high								=	1e12;												% Upper bound for W
    x_best								=	zeros(N, 1);
    tolerance							=	0.1;												% Tolerance for average dB
    max_iter							=	30;													% Maximum bisection iterations
    for iter							=	1:max_iter
        W_mid							=	sqrt(W_low * W_high);

        % Construct real-valued LS system for real x
        sqrtW							=	sqrt(W_mid);
        A_ls							=	[sqrtW * real(Astop); sqrtW * imag(Astop); real(Apass); imag(Apass)];
        b_ls							=	[zeros(2 * size(Astop, 1), 1); ones(size(Apass, 1), 1); zeros(size(Apass, 1), 1)];

        % Solve least squares with constraints using lsqlin
        x								=	lsqlin(A_ls, b_ls, [], [], Aeq, beq, lb, ub, [], optimset('Display', 'off'));

        % Compute average stopband dB
        stopband_mag					=	abs(Astop * x);
        avg_dB							=	mean(20 * log10(stopband_mag + eps));

        if abs(avg_dB - target_avg_dB)	<	tolerance
            x_best						=	x;
            break;
        elseif avg_dB					>	target_avg_dB										% Too high (less rejection), need larger W
            W_low						=	W_mid;
        else																					% Too low, reduce W_high
            W_high						=	W_mid;
            x_best						=	x;
        end
    end

    % Use the best x
    x									=	x_best;
    %x 									=	[1 -2 1];
	
    % Normalize for unity passband gain
    avg_pass							=	mean(Apass * x);
    x									=	x / real(avg_pass);

    % Plot frequency response
    f_plot								=	linspace(0, f_max * 1.1, 5000);
    f_plot_norm							=	f_plot / PRF_avg;
    H									=	zeros(size(f_plot));
    for k								=	1:length(f_plot)
        H(k)							=	x' * exp(-1j * 2 * pi * f_plot_norm(k) * tau_norm');
    end
    figure; hold on; grid on;
    plot(f_plot, 20 * log10(abs(H) + eps));
    plot([0, f_max], [target_avg_dB, target_avg_dB], 'r--');
    xlabel('Doppler Frequency (Hz)'); ylabel('Magnitude (dB)');
    title(['Frequency Response for Sequence Length ' num2str(M)]);
    xlim([0, f_max * 1.1]); ylim([-60, 20]);
    legend('Response', 'Target Average -45 dB');
end