% Step 1: Audio Source 
target_pos = [0, 10];
[target_audio, Fs] = audioread('delightful.wav');

% Step 2: Delay due to distance from audio source

%positions x,y in meters
% center mic, left mic, right mic
mic1 = [0 0];
mic2 = [-1 -1];
mic3 = [1 -1];
mic_pos = [mic1; mic2; mic3]

%distance of each
mic_d = vecnorm((mic_pos - target_pos)');

%delay of each
SOUND = 343; %meters/second
mic_delay = mic_d/SOUND; 

%mic1_d = 

% Model Microphone Transfer Function