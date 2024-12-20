% Constants
c = 343; % Speed of sound (m/s)
f = 1000; % Signal frequency (Hz)
lambda = c / f; % Wavelength (m)
d = 0.01; % Microphone spacing (m)

% Microphone array positions (linear array along x-axis)
M = 4; % Number of microphones
mic_pos = (0:M-1)' * [d, 0, 0]; % Linear array in the x-y plane

% Azimuth angles to compute
azimuth_range = -180:1:180; % Azimuth angles (degrees)
beam_pattern = zeros(1, length(azimuth_range)); % Initialize beam pattern

% Fixed elevation angle
elevation_angle = 0; % Horizontal plane
theta = deg2rad(elevation_angle); % Convert to radians

% Loop over azimuth angles
for az = 1:length(azimuth_range)
    % Compute azimuth unit vector
    phi = deg2rad(azimuth_range(az)); % Convert azimuth to radians
    direction = [cos(theta)*cos(phi), cos(theta)*sin(phi), sin(theta)];
    
    % Compute time delays for each microphone
    delays = (mic_pos * direction') / c;
    
    % Compute steering vector
    steering_vector = exp(-1j * 2 * pi * f * delays);
    
    % Compute beamformer response (dot product of steering vectors)
    beam_pattern(az) = abs(sum(steering_vector));
end

% Normalize beam pattern
beam_pattern = beam_pattern / max(beam_pattern);

% Plot beam pattern
figure;
polarplot(deg2rad(azimuth_range), beam_pattern, 'LineWidth', 1.5);
title('Beam Pattern for Azimuth Angles');
ax = gca;
ax.ThetaZeroLocation = 'top'; % 0° at the top
ax.ThetaDir = 'clockwise'; % Azimuth angles clockwise
grid on;

function y = steering_delay()
    t = (0:10000)/5000;
    X = 0.5*sin(t'*2*pi*262);
    mic_pos = (0:4-1)' * 0.01;
    sv = steervec(mic_pos',[90,0]);
    y = conj(sv) * X;
    disp(sv);
end

%steering_delay();