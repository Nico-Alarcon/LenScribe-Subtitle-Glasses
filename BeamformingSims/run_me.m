clear;clf; close all;
verbose = true;
speak = false;
graph = true;

set(0, 'DefaultFigurePosition', [300, 300, 600, 300]); % Default size for all figures
SOUND = 343; %meters/second

%%  NUMBER AND COORDINATES OF MICROPHONES IN METERS  %

%Endfire Linear Microphone Array
mic_n = 10; % Number of microphones
d = 0.01; % Spacing (1 cm), under nyquist for frequencies of interest 
mic_pos = (0:mic_n-1)' * [0,d]; % Linear positions along y-axis 

theta = pi/2;
target_pos = [cos(theta), sin(theta)];

%% target audio
ns = 60000;
fs_adc = 20000;
ts = (0:ns-1)/fs_adc;
target_audio = 0.5*sin(ts'*2*pi*262);

%% %%%%%% Delay Between Mics %%%%%%
%delta_t, distance/speed
abs_delay = vecnorm((mic_pos - target_pos)')/343;
req_delay = abs_delay;
%   Noise Source  Uris Library%

if graph
    figure;
    subplot(1,2,1);
    %mic positions
    plot(mic_pos(:,1)*10^3, mic_pos(:,2)*10^3, 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    xlim([-60,60]);
    ylim([-60,60]);
    title('Microphone Array Position in mm');
    subplot(1,2,2);
    plot(mic_pos(:,1)*10^3, mic_pos(:,2)*10^3, 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    hold on;
    %target audio
    plot(target_pos(:,1)*10^3, target_pos(:,2)*10^3, 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    hold on;
    % Labels and formatting
    title('Target and Noise Coordinates in mm');
    grid on;
    axis equal; % Ensures equal scaling of x and y axes
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                         SIGNAL CONVERTED TO DIGITAL                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEAMFORMING OPERATIONS
% time domaindelay and sum --> simplest, least computationally intensive
% After computing mic_data and required delay:

signal = zeros(ns, mic_n);
theta = (-180:180)*pi/180;
power_output = zeros(length(theta),1);
for i = 1:length(theta)
    for j = 1:mic_n
        %insert delay # of samples before and attenuate by 4pi*distance
        signal(:,j) = add_source(mic_pos(j,:), [cos(theta(i)), sin(theta(i))], target_audio, fs_adc);
    end
    power_output(i) = bandpower(real(ifft(freq_delaysum(fft(signal), req_delay, fs_adc))), fs_adc, [0,fs_adc/2]);
end
power_output = 10*log10(abs(power_output));
power_output = power_output - max(power_output);


figure;
polarplot(theta, real(power_output));
title('Phase vs. Array Response');
grid on;



%% %%%%%                 HELPER FUNCTIONS                   %%%%%%%%%%%
% function takes source locations, sampling Fs, source locations, and mic position
% outputs resulting audio for that mic
function audio = add_source(input_pos, sources, signal, Fs)
    d = vecnorm((sources - input_pos)'); %distance from mic to noise sources
    delays = int32((d/343)*Fs);
    delays = delays - min(delays);
    audio = zeros(length(signal),1);
    for i = 1:length(delays)
        audio = audio + [zeros(delays(i),1); signal(1:end-delays(i))]/(4*pi*d(i));
    end
    audio = audio/size(sources,1); 
end
