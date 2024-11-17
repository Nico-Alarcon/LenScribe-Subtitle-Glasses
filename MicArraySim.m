clear;

verbose = false;

% Step 1: Audio Source 

target_pos = [0, 2];
[target_audio, Fs] = audioread('recorded_audio.wav');

if verbose
    sound(target_audio, Fs);
end

% Step 2: Delay due to distance from audio source

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

if verbose
    fprintf("The microphones have delays of %s microseconds relative to mic1\n",join(string(abs(diff(mic_delay))*1e6), ', ')); 
end

% Model Microphone Transfer Function, Fs=768kHz 


% Plot Audio in Time and Frequency Domain