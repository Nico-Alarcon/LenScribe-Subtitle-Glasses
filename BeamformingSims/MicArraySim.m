clear;
verbose = true;
speak = false;

%   NUMBER AND POSITIONS OF MICROPHONES IN METERS  %
mic1 = [0 0];
mic2 = [-72.6e-3 -10e-3];
mic3 = [72.6e-3 -10e-3];
mic_pos = [mic1; mic2; mic3];

%   POSITION OF AUDIO SOURCE    %
target_pos = [0, 2];

[target_audio, Fs] = audioread('recorded_audio.wav');
n = length(target_audio);
t = (0:n-1)/Fs;

if speak
    sound(target_audio, Fs);
end

% Step 2: Delay and attenuation due to distance from audio source
% attenuation by 1/(4pi*r)
%positions x,y in meters, rough estimates from my glasses

%distance of each from target audio
mic_d = vecnorm((mic_pos - target_pos)');

% Subplot 1: Time-domain representation
figure;

subplot(4, 1, 1);
plot(t, target_audio);
title('Time-Domain Representation of Audio Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

%delay of each and graph
SOUND = 343; %meters/second
mic_delay = mic_d/SOUND; 
mic_data = zeros(n, length(mic_delay));
for i = 1:length(mic_delay)
    delay = int32(mic_delay(i)*Fs);

    %insert delay # of samples before and attenuate by 4pi*distance
    mic_data(:,i) = [zeros(delay,1); target_audio(1:end-delay)]/(4*pi*mic_d(i));

    subplot(4, 1, i+1);
    plot(t, mic_data(:,i));
    title(sprintf('Mic_{%d}, Delay %f ms', i, round(mic_delay(i)*1e3,3)));
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
end

if verbose
    fprintf("The microphones have delays of %s microseconds relative to mic1\n",join(string(abs(diff(mic_delay(1:2)))*1e6), ', ')); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       SIGNAL REACHES MICROPHONES                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MICROPHONE CONTINUOUS TIME H(s)


% ANTIALIASING FILTER


% SAMPLING


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         SIGNAL CONVERTED TO DIGITAL                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BEAMFORMING OPERATIONS

% delay and sum --> simplest, least computationally intensive

% Bartlett Beamforming --> delay/sum + weighting per microphone

% frequency domain delay/sum + deconvolution --> more computation, cleaner


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         COMPARISONS AND ANALYSIS                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% COMPARISON BETWEEN EACH


% fft audio data, multiply amplitude add phase

% Frequency points (Hz)
f = [20 40 60 80 100 200 400 600 800 1000 2000 4000 6000 8000 10000 20000];
% Amplitude response in dB (relative to a reference, e.g., 1kHz)
magnitude_dB = [-3 -1.25 -0.8 -0.40 -0.39 -0.05 0 0 0 0 0 0.01 0.18 0.2 1 4];
% Phase response (degrees)
phase_deg = [62 37 22 18 16 8 5 3 1 0.5 0 -1 -2 -4.8 -5.1 -5.2];

% Call the function
[magnitudeCoeffs, phaseCoeffs] = fitFrequencyResponse(f, magnitude_dB, phase_deg, 5);
evalFunction(magnitudeCoeffs, 20)





%%%%%%%                 HELPER FUNCTIONS                   %%%%%%%%%%%
function out = evalFunction(coeff, f)
    %from an array of coefficients evaluate function with log input
    out = 0;
    for i = 1:length(coeff)
        out = out + coeff(i)*log10(f)^(length(coeff)-i);
    end
end