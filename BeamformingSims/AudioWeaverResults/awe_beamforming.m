%% Four-Microphone Array Generator in Millimeters with CSV Export, 0° Arrow, and Formatted Console Output
% Three microphones are positioned on a circle of radius 0.04 m (40 mm)
% at angles 0°, 120°, and 240°. These positions are fixed.
%
% The fourth microphone ("origin mic") is normally at the origin,
% but here its x-coordinate is tunable. Its default value is set to half of the
% 120° mic’s x-coordinate. For a 0.04 m radius:
%   0.04 * cosd(120) = 0.04 * (-0.5) = -0.02 m, so half of that is -0.01 m.
%
% After computing the positions in meters, the coordinates are converted
% to millimeters. They are then exported in two formats:
%   1. "mic_array_doa_sim.txt" contains the coordinates as rows ([x, y] per row).
%   2. "mic_array_find_doa.txt" contains the coordinates transposed so that
%      columns contain the x and y coordinates.
%
% Additionally, the script prints the x- and y-coordinate arrays in the format:
%   [x1, x2, x3, ...]
%   [y1, y2, y3, ...]
% (All values are in millimeters.)

%% Parameters
circleRadius = 0.015;  % Default circle radius in meters

% Calculate the 120° mic's x-coordinate
mic120_x = circleRadius * cosd(120);  % For 120°, cosd(120) = -0.5 => mic120_x = -0.02 m

% Default origin mic x-coordinate: half of the 120° mic's x-coordinate.
defaultOriginMicX = mic120_x / 2;      % Defaults to -0.01 m
originMicX = defaultOriginMicX;        % Tunable origin mic x-coordinate (default -0.01 m)

%% Compute microphone positions in meters

% Microphones on the circle (fixed positions)
mic0   = [circleRadius, 0];                           % 0° mic: (0.04, 0)
mic120 = [circleRadius * cosd(120), circleRadius * sind(120)];  % 120° mic
mic240 = [circleRadius * cosd(240), circleRadius * sind(240)];  % 240° mic

% Origin mic with tunable x (y remains 0)
originMic = [originMicX, 0];

% Combine positions into a matrix (each row is [x, y] in meters)
micPositions = [mic0; mic120; mic240; originMic];

%% Convert positions to millimeters
micPositions_mm = micPositions * 1000;  % Conversion: 1 m = 1000 mm

%% Display the positions in millimeters
disp('Generated Microphone Positions (in millimeters):');
fprintf('0° Mic:     [%.2f, %.2f]\n', micPositions_mm(1,1), micPositions_mm(1,2));
fprintf('120° Mic:   [%.2f, %.2f]\n', micPositions_mm(2,1), micPositions_mm(2,2));
fprintf('240° Mic:   [%.2f, %.2f]\n', micPositions_mm(3,1), micPositions_mm(3,2));
fprintf('Origin Mic: [%.2f, %.2f] (Default x = half of 120° mic x)\n', micPositions_mm(4,1), micPositions_mm(4,2));

%% Extract separate arrays for x and y coordinates (in mm)
x_coords = micPositions_mm(:,1);
y_coords = micPositions_mm(:,2);

% Format the coordinate arrays as strings:
x_str = ['[', sprintf('%.2f, ', x_coords)];
x_str = [x_str(1:end-2), ']'];  % Remove trailing comma/space, add closing bracket

y_str = ['[', sprintf('%.2f, ', y_coords)];
y_str = [y_str(1:end-2), ']'];

% Print the formatted coordinate arrays to the console
disp('X Coordinates (mm):');
disp(x_str);
disp('Y Coordinates (mm):');
disp(y_str);

%% Define file paths (adjust if necessary)
basePath = 'C:\DSP Concepts\AWE Designer 8.D.2.7 Standard';
filePathRows = fullfile(basePath, 'mic_array_doa_sim.txt');
filePathCols = fullfile(basePath, 'mic_array_find_doa.txt');

%% Export the coordinates (in mm)
% 1. Write the micPositions_mm matrix as rows ([x, y] per line)
writematrix(micPositions_mm, filePathRows);
fprintf('Microphone coordinates exported (rows) to %s\n', filePathRows);

% 2. Write the transposed matrix so that columns contain x and y coordinates.
%    This creates a 2x4 matrix: first row contains x values, second row contains y values.
writematrix(micPositions_mm.', filePathCols);
fprintf('Microphone coordinates exported (columns) to %s\n', filePathCols);

%% Plot the microphone array for visualization (in mm)
figure;
hold on;
% Plot the three circle mics
plot(micPositions_mm(1,1), micPositions_mm(1,2), 'bo', 'MarkerSize', 10, 'LineWidth', 2);
plot(micPositions_mm(2,1), micPositions_mm(2,2), 'bo', 'MarkerSize', 10, 'LineWidth', 2);
plot(micPositions_mm(3,1), micPositions_mm(3,2), 'bo', 'MarkerSize', 10, 'LineWidth', 2);
% Plot the origin mic (tunable x)
plot(micPositions_mm(4,1), micPositions_mm(4,2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
% Draw the reference circle (converted to mm)
theta = linspace(0, 2*pi, 200);
x_circle = (circleRadius * cos(theta)) * 1000;
y_circle = (circleRadius * sin(theta)) * 1000;
plot(x_circle, y_circle, 'k--');

xlabel('X (mm)');
ylabel('Y (mm)');
title('Four-Microphone Array Configuration (mm)');
legend('0° Mic', '120° Mic', '240° Mic', 'Origin Mic (tunable)', 'Reference Circle');
axis equal;
grid on;
hold off;

%% Add an arrow annotation for 0° direction
% The annotation uses normalized figure coordinates (values between 0 and 1).
annotation('textarrow',[0.3, 0.6],[0.5, 0.5],'String','0 degrees','FontSize',12);
