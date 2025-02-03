% Clear workspace and figures
clear; close all; clc;

% Parameters
c = 343;                % Speed of sound (m/s)
fs = 44100;             % Sampling frequency (Hz)
t = 0:1/fs:0.1;         % Time vector (0.1 seconds)
f0 = 2000;              % Source signal frequency (Hz)
noise_level = 0.2;      % Noise level (standard deviation)

% Microphone array geometry (tetrahedron)
mic_positions = [1, 1, 1;           % Vertex 1
                -1, -1, 1;          % Vertex 2
                -1, 1, -1;          % Vertex 3
                1, -1, -1];         % Vertex 4
mic_positions = mic_positions * 0.05; % Scale to 5cm edge length

% Source direction (azimuth and elevation in degrees)
azimuth = 45;
elevation = 30;

% Convert spherical coordinates to Cartesian unit vector
az = deg2rad(azimuth);
el = deg2rad(elevation);
u = [cos(el)*cos(az); cos(el)*sin(az); sin(el)]; % DOA unit vector

% Calculate time delays for each microphone
tau = mic_positions * u / c; % Propagation delays (seconds)

% Generate source signal (sine wave)
s = sin(2*pi*f0*t);

% Simulate received signals with delays and noise
received_signals = zeros(4, length(t));
for i = 1:4
    % Apply time delay using frequency-domain phase shift
    delay_phase = exp(-1i*2*pi*f0*tau(i));
    delayed_s = real(s * delay_phase);
    
    % Add Gaussian noise
    received_signals(i,:) = delayed_s + noise_level*randn(size(delayed_s));
end

% Frequency-domain beamforming
nfft = 2^nextpow2(length(t));       % FFT length
frequencies = (0:nfft-1)*(fs/nfft); % Frequency bins

% Calculate steering vector
steering_vector = exp(-1i*2*pi*frequencies'*tau');

% FFT of received signals
S = fft(received_signals, nfft, 2);

% Apply beamforming weights and sum
beamformed_spectrum = sum(S .* conj(steering_vector), 1);

% Inverse FFT to get time-domain signal
beamformed_signal = ifft(beamformed_spectrum, 'symmetric');
beamformed_signal = beamformed_signal(1:length(t));

% Calculate SNR improvement
original_snr = 20*log10(norm(s)/norm(s - received_signals(1,:)));
beamformed_snr = 20*log10(norm(s)/norm(s - beamformed_signal));

% Plot results
figure;
subplot(3,1,1);
plot(t, s);
title('Original Source Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(3,1,2);
plot(t, received_signals(1,:));
title(['Received Signal at Mic 1 (SNR = ' num2str(round(original_snr,1)) ' dB)']);
xlabel('Time (s)');
ylabel('Amplitude');

subplot(3,1,3);
plot(t, beamformed_signal);
title(['Beamformed Signal (SNR = ' num2str(round(beamformed_snr,1)) ' dB)']);
xlabel('Time (s)');
ylabel('Amplitude');

% Plot array geometry
figure;
scatter3(mic_positions(:,1), mic_positions(:,2), mic_positions(:,3), 100, 'filled');
hold on;
quiver3(0, 0, 0, u(1), u(2), u(3), 0.2, 'r', 'LineWidth', 2);
title('4-Mic Trillium Array Geometry');
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
legend('Microphones', 'DOA Direction');
grid on; axis equal;