clear;clf; close all;
verbose = false;
speak = false;

set(0, 'DefaultFigurePosition', [100, 100, 1200, 800]); % Default size for all figures

%   NUMBER AND POSITIONS OF MICROPHONES IN METERS  %

%Endfire Linear Microphone Array
mic_n = 4; % Number of microphones
d = 0.01; % Spacing (1 cm), under nyquist for frequencies of interest 
mic_pos = (0:mic_n-1)' * [0,d]; % Linear positions along y-axis 

%   POSITION OF AUDIO SOURCE    %
target_pos = [0, 2];
[target_audio, Fs] = audioread('recorded_audio.wav');
n = length(target_audio);
t = (0:n-1)/Fs;
f = (0:n-1)*(Fs/n); % Frequency axis (Hz)

% pure tone overwrite for testing, middle C
target_audio = 0.5*sin(t'*2*pi*262);

%   Noise Source  %
[noise_audio, Fs2] = audioread('AmbientNoise.wav');

%resample if different
[p,q] = rat(Fs/Fs2);
noise_audio = resample(noise_audio, p, q);
noise_audio = noise_audio(1:n)*2; %make a bit louder

%% positions of noise source, roughly parabolic with noise 2 feet behind and
% to sides 2x^2 - 2, mics at 0,0
noise_x = 2; % length is noise_x*2 + 1
noise_pos = horzcat((-noise_x:noise_x)',((-noise_x:noise_x)'.^2).*2 -2);

if speak
    disp("Target Audio");
    sound(target_audio, Fs);
    pause(length(target_audio)/Fs +1); %stop for another second
end

%% Step 2: Delay and attenuation due to distance from audio source
% attenuation by 1/(4pi*r) as sound pressure is the metric of interest


%% Subplot 1: Time-domain representation
figure;

subplot(mic_n+1, 1, 1);
plot(t, target_audio, 'DisplayName', 'Target', 'LineWidth', 1.5); hold on;
plot(t, noise_audio,  'DisplayName', 'Noise', 'LineWidth', 1.5);

title('Time-Domain: Point Source Target vs. Noise');
legend('Location', 'best');
xlabel('Time (s)');
xlim([0,max(t)])
ylabel('Amplitude');
grid on;

%%delay of each and graph
SOUND = 343; %meters/second
mic_d = vecnorm((mic_pos - target_pos)'); %distance
mic_delay = mic_d/SOUND; %delay in time
mic_data = zeros(n, mic_n);
signal = zeros(n, mic_n);
for i = 1:mic_n
    %insert delay # of samples before and attenuate by 4pi*distance
    signal(:,i) = add_source(mic_pos(i,:), target_pos, target_audio, Fs);
    mic_data(:,i) = signal(:,i) + add_source(mic_pos(i,:), noise_pos, noise_audio, Fs);

    subplot(mic_n+1, 1, i+1);
    plot(t, mic_data(:,i));
    title(sprintf('Mic_{%d}, Delay %f ms', i, round(mic_delay(i)*1e3,3)));
    xlim([0,max(t)])
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on; hold on;
end

% target audio + noise as it reaches microphones
if speak 
    disp("Mic Input");
    sound(mean(mic_data,2),Fs);
    pause(length(mic_data(:,1))/Fs + 1);
end

mic_rel_delay = abs(mic_delay - max(mic_delay));
if verbose
    for i = 1:length(mic_delay)
        fprintf("Mic %d needs a delay of %sus\n",i, string(mic_rel_delay(i)*1e6)); 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                       SIGNAL REACHES MICROPHONES                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% MICROPHONE CONTINUOUS TIME |H(s)|  
f_approx = [35, 100, 900, 1000, 1100, 8000, 15000]; %Hz
magnitude_dB = [-1.0, -0.5, 0.0, 0.0, 0.0, 1.5, 4.5];
phase_deg = [0, 0, 0, 0, 0, 0, 0];
% Mic Frequency Response
mic_H = fitFrequencyResponse(f_approx, magnitude_dB, phase_deg, 5, f);

% ANTIALIASING FILTER, 4th order bessel filter to preserve phase
N=4;
fp = 2*pi*8000; %cutoff frequency
[z, p, k] = besself(N,fp);
[num, den] = zp2tf(z, p, k);
[num_d, den_d] = impinvar(num, den, Fs); %convert to digital
[H,f] = freqz(num_d,den_d, n, Fs);

%% Plotting Group Delay and Magnitude Response of Filter
figure;

subplot(2,1,1);
semilogx(f,20*log(abs(H)));
title("Magnitude Response of AntiAliasing Filter")
xlabel("Frequency (Hz)")
ylabel("Magnitude (dB)")
xline(fp)
grid on

subplot(2,1,2);
grpdel = -diff(unwrap(angle(H)))./diff(2*pi*f);
semilogx(f(2:end),grpdel)
title("Group Delay of AntiAliasing Filter")
xlabel("Frequency (Hz)")
ylabel("Group delay (s)")
xline(fp)
grid on


%% %%%%%%%%% Apply Microphone and AAF tfs %%%%%%%%%%%%%%%%%%
input_H = H .* mic_H;
beam_input = fft(mic_data) .* repmat(input_H,1,mic_n);

%%%%%%%%%%% FFT Graphs  %%%%%%%%%%%%%%%% Should Look V Similar
audio_fft = abs(fft(target_audio)); % Compute the FFT

figure;
subplot(2,1,1);
semilogx(f(1:n/2), 20*log10(audio_fft(1:n/2,1))); % Plot only the positive frequencies
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
xlim([f(1),f(n/2)]);
ylim([-200, 100]);
title('Magnitude Spectrum of Target Audio Signal');
subplot(2,1,2);

semilogx(f(1:n/2), 20*log10(abs(beam_input(1:n/2))), 'DisplayName','Beamforming Input'); hold on;
semilogx(f(1:n/2), 20*log10(abs(H(1:n/2))), 'b', 'LineWidth', 2,'DisplayName', 'AAF Filter');
semilogx(f(1:n/2), 20*log10(abs(mic_H(1:n/2))), 'r' , 'LineWidth', 2, 'DisplayName', 'Microphone H(s)');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
legend('Location', 'best');
xlim([f(1),f(n/2)]);
ylim([-200, 100]);
title('Magnitude Spectrum of Audio Signal After Mics and AAF');

%% %%%%%%%%%%%%% FFT Check %%%%%%%%%%%%%

% SAMPLING - uses its own AAF, but doesnt matter if fs > 16kHz
fs_adc = 20000;
[p, q] = rat(fs_adc / Fs); % Calculate resampling factors

beam_input = real(ifft(beam_input));
beam_input = resample(beam_input, p, q); % Resample signal
ns = length(beam_input);
ts = (0:ns-1)/fs_adc;
f_resample = linspace(0, fs_adc, ns);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                         SIGNAL CONVERTED TO DIGITAL                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEAMFORMING OPERATIONS
phi = 90*pi/180;
req_delay = (mic_pos * [cos(phi),sin(phi)]'/SOUND)'; %same as mic_rel_delay

%% delay and sum --> simplest, least computationally intensive
 % After computing mic_data and mic_delay:
y_beamformed = delaysum(beam_input, req_delay, fs_adc);
audiowrite('delaysum_signal.wav', y_beamformed, fs_adc);

%% Bartlett Beamforming = delay and suum in frequency domain
%same delay for same distance and desired endfire response
delaysum_bartlett_fft = bartlett(fft(beam_input), f_resample, req_delay);
delaysum_steer_time = real(ifft(delaysum_bartlett_fft));

if speak
    disp("Delay and Sum Beamforming");
    sound(y_beamformed, fs_adc);
    pause(length(y_beamformed)/Fs +1);

    disp("Frequency Domain Delay and Sum Beamforming");
    sound(delaysum_steer_time, fs_adc);
    pause(length(delaysum_steer_time)/fs_adc)
end

% Plotting for Comparison 
% Plot time-domain: Overlay beamformed signal with the original
figure;
subplot(2,1,1);

plot(ts, y_beamformed, 'DisplayName', 'Beamformed', 'LineWidth', 1.5);
hold on;
plot(ts, delaysum_steer_time, 'DisplayName', 'Freq Domain Beamformed', 'LineWidth', 1.5);
hold on;
plot(ts, mean(beam_input,2), 'DisplayName', 'Pre-Beamformed', 'LineWidth', 1.5);
hold on;
plot(t, mean(signal,2), 'DisplayName', 'Original', 'LineWidth', 1.5);
hold off;

title('Time-Domain: Original vs Beamformed');
xlabel('Time (s)');
ylabel('Amplitude');
xlim([t(1),t(end)]);
legend('Location', 'best');
grid on;

% Frequency-domain comparison
Y_original = fft(mean(signal,2));
Y_prebeamformed = fft(mean(beam_input,2));
Y_beamformed = fft(y_beamformed);

% Convert to magnitude (dB) for better comparison
Y_original_mag_dB = 20*log10(abs(Y_original(1:n/2)));
Y_beamformed_mag_dB = 20*log10(abs(Y_beamformed(1:ns/2)));
Y_prebeamformed_mag_dB = 20*log10(abs(Y_prebeamformed(1:ns/2)));
Y_delaysum_steer_mag_dB = 20*log10(abs(delaysum_bartlett_fft(1:ns/2)));

subplot(2,1,2);
semilogx(f_resample(1:ns/2), Y_beamformed_mag_dB, 'DisplayName', 'Beamformed', 'LineWidth', 1.5);
hold on;
semilogx(f_resample(1:ns/2), Y_prebeamformed_mag_dB, 'DisplayName', 'Pre-Beamformed', 'LineWidth', 1.5);
hold on;
semilogx(f(1:n/2), Y_original_mag_dB, 'DisplayName', 'Original', 'LineWidth', 1.5);
hold on;
semilogx(f_resample(1:ns/2), Y_delaysum_steer_mag_dB , 'DisplayName', 'Freq Domain Beamforming', 'LineWidth', 1.5);
hold off;

title('Frequency-Domain: Original vs Beamformed');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
legend('Location', 'best');

xlim([f(1),min([f(ns/2),f_resample(ns/2)])]);
ylim([-200, 100]);
grid on;


%plot(f_resample, delaysum_steer_fft, 'DisplayName', 'Freq Domain Beamforming', 'LineWidth', 1.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                         COMPARISONS AND ANALYSIS                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% COMPARISON BETWEEN EACH


% fft audio data, multiply amplitude add phase

% Frequency points (Hz)




%%%%%%%                 HELPER FUNCTIONS                   %%%%%%%%%%%
% function takes source signal, sampling Fs, source locations, and mic position
% outputs resulting audio
function audio = add_source(input_pos, sources, signal, Fs)
    d = vecnorm((sources - input_pos)'); %distance from mic to noise sources
    delays = int32((d/343)*Fs); 
    audio = zeros(length(signal),1);
    for i = 1:length(delays)
        audio = audio + [zeros(delays(i),1); signal(1:end-delays(i))]/(4*pi*d(i));
    end
    audio = audio/size(sources,1); 
end

