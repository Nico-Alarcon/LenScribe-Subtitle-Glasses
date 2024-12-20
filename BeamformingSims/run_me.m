clear;clf; close all;
verbose = true;
speak = true;
graph = true;

set(0, 'DefaultFigurePosition', [300, 300, 600, 300]); % Default size for all figures
SOUND = 343; %meters/second

%%  NUMBER AND COORDINATES OF MICROPHONES IN METERS  %

%Endfire Linear Microphone Array
mic_n = 4; % Number of microphones
d = 0.01; % Spacing (1 cm), under nyquist for frequencies of interest 
mic_pos = (0:mic_n-1)' * [0,d]; % Linear positions along y-axis 
target_pos = [0, 2];

%   POSITION OF AUDIO SOURCE   3 seconds of audio %
Fs = 44100; %from testing with audiorecorder
n = 3*Fs;
f = linspace(0,Fs,n); % Frequency axis (Hz)
t = (0:n-1)/Fs;
% pure tone overwrite for testing, middle C
target_audio = 0.5*sin(t'*2*pi*262);


%%%%%%%% Delay Between Mics %%%%%%
phi = 90*pi/180;
%same relative delay between them = distance between mics and speed constant
req_delay = (mic_pos * [cos(phi),sin(phi)]'/SOUND)';

abs_delay = vecnorm((mic_pos - target_pos)')/343;
%   Noise Source  Uris Library%
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

% Figure 1: Time-domain representation of Signals at Mics
if graph
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
end
%%delay of each and graph



mic_data = zeros(n, mic_n);
signal = zeros(n, mic_n);
for i = 1:mic_n
    %insert delay # of samples before and attenuate by 4pi*distance
    signal(:,i) = add_source(mic_pos(i,:), target_pos, target_audio, Fs);
    mic_data(:,i) = signal(:,i) + add_source(mic_pos(i,:), noise_pos, noise_audio, Fs);
    
    if graph
        subplot(mic_n+1, 1, i+1);
        plot(t, mic_data(:,i));
        title(sprintf('Mic_{%d}, Delay %f ms', i, round(abs_delay(i)*1e3,3)));
        xlim([0,max(t)])
        xlabel('Time (s)');
        ylabel('Amplitude');
        grid on; hold on; 
    end
end

% target audio + noise as it reaches microphones
if speak 
    disp("Mic Input");
    sound(mean(mic_data,2),Fs);
    pause(length(mic_data(:,1))/Fs + 1);
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
H = freqz(num_d,den_d, f, Fs);


%% Plotting Group Delay and Magnitude Response of Filter
% Figure 2: Bessel AntiAliasing Filter, Cutoff at 8kHz
if graph
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
end

%% %%%%%%%%% Apply Microphone and AAF tfs %%%%%%%%%%%%%%%%%%
input_H = H.' .* mic_H;
beam_input = fft(mic_data) .* repmat(input_H,1,mic_n);

%%%%%%%%%%% FFT Graphs  %%%%%%%%%%%%%%%% Should Look V Similar
audio_fft = abs(fft(target_audio)); % Compute the FFT

%Figure 3: Target Audio Signal FFT and Target + Noise
if graph
    figure;
    subplot(2,1,1);
    semilogx(f(1:n/2), 20*log10(audio_fft(1:n/2,1))); % Plot only the positive frequencies
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    xlim([f(1),f(n/2)]);
    ylim([-200, 100]);
    title('Magnitude Spectrum of Target Audio Signal Without Noise');
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
end

%% %%%%%%%%%%%%% FFT Check %%%%%%%%%%%%%
% SAMPLING - C/D Converstion - uses its own AAF, but we did our own first
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
beam_input_fft = fft(beam_input);
% time domaindelay and sum --> simplest, least computationally intensive
% After computing mic_data and required delay:
delaysum_time_time = time_delaysum(beam_input, req_delay, fs_adc);
audiowrite('delaysum_signal.wav', delaysum_time_time, fs_adc);

%% delay and sum in frequency domain
%same delay for same distance and desired endfire response
delaysum_freq_fft = freq_delaysum(beam_input_fft, f_resample, req_delay);
delaysum_freq_time = real(ifft(delaysum_freq_fft));

%% Bartlett Beamforming
theta = (-90:90)*pi/180;
power_output = zeros(length(theta),1);
for i = 1:length(theta)
    [power_output(i), ~] = bartlett(beam_input, d, mic_n, theta(i), fs_adc);
end
power_output = 10*log10(abs(power_output));
power_output = power_output - max(power_output);

[~,y_bartlett] = bartlett(beam_input, d, mic_n, 0*pi/180, fs_adc);
y_bartlett = real(y_bartlett);
% :(
if graph
    figure;
    plot(theta*180/pi, real(power_output));
    title('Phase vs. Bartlett Array Response');
    xlabel('Phase (degrees)');
    ylabel('Normalized Power Output (dB)');
    xlim([-90,90]);
    grid on;
end

if speak
    disp("Delay and Sum Beamforming");
    sound(delaysum_time_time, fs_adc);
    pause(length(delaysum_time_time)/Fs +1);

    disp("Frequency Domain Delay and Sum Beamforming");
    sound(delaysum_freq_time, fs_adc);
    pause(length(delaysum_freq_time)/fs_adc);

    disp("Bartlett Beamforming");
    sound(y_bartlett,fs_adc);
    pause(length(y_bartlett)/fs_adc);
end

if verbose
    for i = 1:length(req_delay)
        fprintf("Mic %d needs a delay of %sus\n",i, string(req_delay(i)*1e6)); 
    end
end

% Frequency-domain comparison
Y_original = fft(target_audio);
Y_prebeamformed = fft(beam_input(:,1));
Y_delaysum_time_fft = fft(delaysum_time_time);
Y_bartlett = fft(y_bartlett);

% Convert to magnitude (dB) for better comparison
Y_original_mag_dB = 20*log10(abs(Y_original(1:n/2)));
Y_delaysum_time_mag_dB = 20*log10(abs(Y_delaysum_time_fft(1:ns/2)));
Y_prebeamformed_mag_dB = 20*log10(abs(Y_prebeamformed(1:ns/2)));
Y_delaysum_freq_mag_dB = 20*log10(abs(delaysum_freq_fft(1:ns/2)));
Y_bartlett_dB = 20*log10(abs(Y_bartlett(1:ns/2)));

% Plotting for Comparison 
% Plot time-domain: Overlay beamformed signal with the original
if graph    
    figure;
    subplot(2,1,1);
    
    plot(ts, delaysum_time_time, 'DisplayName', 'Time Delay and Sum', 'LineWidth', 1.5);
    hold on;
    plot(ts, delaysum_freq_time, 'DisplayName', 'Freq Delay and Sum', 'LineWidth', 1.5);
    hold on;
    plot(ts, y_bartlett,'DisplayName', 'Bartlett Beamformed', 'LineWidth', 1.5);
    hold on;
    plot(ts, beam_input(:,1), 'DisplayName', 'Pre-Beamformed', 'LineWidth', 1.5);
    hold on;
    plot(t, target_audio, 'DisplayName', 'Original', 'LineWidth', 1.5);
    hold off;
    
    title('Time-Domain: Original vs Beamformed');
    xlabel('Time (s)');
    ylabel('Amplitude');
    xlim([t(1),t(end)]);
    legend('Location', 'best');
    grid on;
    
    
    subplot(2,1,2);
    semilogx(f_resample(1:ns/2), Y_delaysum_time_mag_dB, 'DisplayName', 'Time Delay and Sum', 'LineWidth', 1.5);
    hold on;
    semilogx(f_resample(1:ns/2), Y_delaysum_freq_mag_dB , 'DisplayName', 'Freq Delay and Sum', 'LineWidth', 1.5);
    hold on;
    semilogx(f_resample(1:ns/2), Y_bartlett_dB , 'DisplayName', 'Bartlett Beamformed', 'LineWidth', 1.5);
    hold on;
    semilogx(f_resample(1:ns/2), Y_prebeamformed_mag_dB, 'DisplayName', 'Pre-Beamformed', 'LineWidth', 1.5);
    hold on;
    semilogx(f(1:n/2), Y_original_mag_dB, 'DisplayName', 'Original', 'LineWidth', 1.5);
    hold off;

    title('Frequency-Domain: Original vs Beamformed');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    legend('Location', 'best');
    
    xlim([f(1),min([f(ns/2),f_resample(ns/2)])]);
    ylim([-200, 100]);
    grid on;

    [~,y_bartlett2] = bartlett(beam_input, d, mic_n, 90*pi/180, fs_adc);
    y_bartlett2 = real(y_bartlett2);
    figure;
    subplot(2,1,1);
    plot(ts, y_bartlett2,'DisplayName', 'Bartlett 90 degrees', 'LineWidth', 1.5);
    hold on;
    plot(ts, y_bartlett,'DisplayName', 'Bartlett 0 degrees', 'LineWidth', 1.5);
    hold on;
    title('Time-Domain: Bartlett 0 deg vs Bartlett 90 deg');
    xlabel('Time (s)');
    ylabel('Amplitude');
    xlim([t(1),t(end)]);
    legend('Location', 'best');
    grid on;

    subplot(2,1,2);
    semilogx(f_resample(1:ns/2), Y_bartlett_dB , 'DisplayName', 'Bartlett 0 degrees', 'LineWidth', 1.5);
    hold on;
    Y_bartlett2 = fft(y_bartlett2);
    Y_bartlett2_dB = 20*log10(abs(Y_bartlett2(1:ns/2)));
    semilogx(f_resample(1:ns/2), Y_bartlett2_dB , 'DisplayName', 'Bartlett 90 degrees', 'LineWidth', 1.5);
    hold on;

    title('Frequency-Domain: Bartlett 0 deg vs Bartlett 90 deg');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    legend('Location', 'best');
    
    xlim([f(1),min([f(ns/2),f_resample(ns/2)])]);
    ylim([-200, 100]);
    grid on;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                         COMPARISONS AND ANALYSIS                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPARISON BETWEEN EACH
%need to resample target audio for comparison
target_var_temp = fft(noise_audio + target_audio);
target_var_temp = abs(target_var_temp/max(target_var_temp)).^2;

[p,q] = rat(fs_adc/Fs);
target_audio = resample(target_audio, p, q);

time_var_temp = abs(Y_delaysum_time_fft/max(Y_delaysum_time_fft)).^2;
freq_var_temp = abs(delaysum_freq_fft/max(delaysum_freq_fft)).^2;
bart_var_temp = abs((Y_bartlett/max(Y_bartlett))).^2;

%% 
if verbose
    fprintf('\nVariance at Input: %s\n', var(target_var_temp(1:n/2)));
    fprintf('Variance at Time Delay Sum Output: %s\n', var(time_var_temp(1:ns/2)));
    fprintf('Variance at Freq Delay Sum Output: %s\n', var(freq_var_temp(1:ns/2)));
    fprintf('Variance at Bartlett Output: %s\n', var(bart_var_temp(1:ns/2)));
end
%% %%%%%                 HELPER FUNCTIONS                   %%%%%%%%%%%
% function takes source locations, sampling Fs, source locations, and mic position
% outputs resulting audio for that mic
function audio = add_source(input_pos, sources, signal, Fs)
    d = vecnorm((sources - input_pos)'); %distance from mic to noise sources
    delays = int32((d/343)*Fs); 
    audio = zeros(length(signal),1);
    for i = 1:length(delays)
        audio = audio + [zeros(delays(i),1); signal(1:end-delays(i))]/(4*pi*d(i));
    end
    audio = audio/size(sources,1); 
end
