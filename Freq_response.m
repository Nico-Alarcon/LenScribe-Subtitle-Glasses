% Frequency points (Hz)
f = [20 40 60 80 100 200 400 600 800 1000 2000 4000 6000 8000 10000 20000]; % Example points; refine based on actual data
% Amplitude response in dB (relative to a reference, e.g., 1kHz)
magnitude_dB = [-3 -1.25 0 0 0 1 4]; % Example values, use actual response data
% Phase response (degrees)
phase_deg = [62 16 0 -2 -5 -6]; % Example values; use actual response data

% Convert to linear scale
magnitude = 10.^(magnitude_dB / 20);

% Convert phase to radians
phase_rad = deg2rad(phase_deg);

% Frequency to angular frequency
omega = 2 * pi * f;

% Create the complex transfer function
H = magnitude .* exp(1j * phase_rad);

% Plot Frequency Response (Optional)
figure;
subplot(2,1,1);
semilogx(f, magnitude_dB, 'b-o'); % Magnitude in dB
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;
title('Magnitude Response');

subplot(2,1,2);
semilogx(f, phase_deg, 'r-o'); % Phase in degrees
xlabel('Frequency (Hz)');
ylabel('Phase (degrees)');
grid on;
title('Phase Response');

% Convert to time domain or define a transfer function
% Example: Interpolate and create a continuous system
f_interp = logspace(log10(20), log10(20000), 500); % Interpolated frequency range
H_interp = interp1(f, H, f_interp, 'spline'); % Interpolated response

% Define a transfer function (optional approximation)
[b, a] = invfreqz(H_interp, 2 * pi * f_interp, 10, 10); % Approximate numerator and denominator
sys = tf(b, a);

% Frequency Response of the approximated system
figure;
freqz(b, a, 1024, 44100); % Example with 44.1 kHz sampling rate
title('Approximated Transfer Function');
