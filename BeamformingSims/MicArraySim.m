clear;

verbose = true;
speak = false;
% Step 1: Audio Source 

target_pos = [0, 1];
[target_audio, Fs] = audioread('recorded_audio.wav');

if speak
    sound(target_audio, Fs);
end

% Step 2: Delay and attenuation due to distance from audio source
% attenuation by 1/(4pi*r)
%positions x,y in meters, rough estimates from my glasses
% center mic, left mic, right mic
mic1 = [0 0];
mic2 = [-72.6e-3 -10e-3];
mic3 = [72.6e-3 -10e-3];
mic_pos = [mic1; mic2; mic3];

%distance of each
mic_d = vecnorm((mic_pos - target_pos)');

%delay of each
SOUND = 343; %meters/second
mic_delay = mic_d/SOUND; 

%sample delay: seconds * Fs
mic_data = mic_delay * Fs;
if verbose
    fprintf("The microphones have delays of %s microseconds relative to mic1\n",join(string(abs(diff(mic_delay(1:2)))*1e6), ', ')); 
end

% Model Microphone Transfer Function, Fs=768kHz 
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


% Plot Audio in Time and Frequency Domain
% Step 1: Load the audio file

target_audio = target_audio(:, 1);  % Use only the first channel if stereo

% Step 2: Define time vector
n = length(target_audio);     % Number of samples
t = (0:n-1) / Fs;       % Time vector

% Step 3: Compute the FFT of the signal
y_fft = fft(target_audio);    % Perform FFT
P2 = abs(y_fft/n);      % Two-sided spectrum
P1 = P2(1:n/2+1);       % One-sided spectrum
P1(2:end-1) = 2*P1(2:end-1); % Adjust amplitude for one-sided spectrum

% Step 4: Define the frequency axis
f = Fs * (0:(n/2)) / n; % Frequency vector

% Step 5: Plot the time-domain signal and its frequency spectrum
figure;

% Subplot 1: Time-domain representation
subplot(2, 1, 1);
plot(t, target_audio);
title('Time-Domain Representation of Audio Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Subplot 2: Frequency-domain representation
subplot(2, 1, 2);
plot(f, P1);
title('Frequency-Domain Representation (FFT) of Audio Signal');
xlabel('Frequency (Hz)');
ylabel('|P1(f)|');
grid on;



%%%%%%%                 HELPER FUNCTIONS                   %%%%%%%%%%%
function out = evalFunction(coeff, f)
    %from an array of coefficients evaluate function with log input
    out = 0;
    for i = 1:length(coeff)
        out = out + coeff(i)*log10(f)^(length(coeff)-i);
    end
end