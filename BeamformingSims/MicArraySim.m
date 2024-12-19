clear;clf;
verbose = false;
speak = false;

%   NUMBER AND POSITIONS OF MICROPHONES IN METERS  %

%Endfire Linear Microphone Array
mic_n = 4; % Number of microphones
d = 0.01; % Spacing (1 cm), under nyquist for frequencies of interest 
mic_pos = horzcat(zeros(mic_n,1),(0:mic_n-1)' * d); % Linear positions along x-axis

%   POSITION OF AUDIO SOURCE    %
target_pos = [0, 2];
[target_audio, Fs] = audioread('recorded_audio.wav');
n = length(target_audio);
t = (0:n-1)/Fs;
f = (0:n-1)*(Fs/n); % Frequency axis (Hz)

%   Noise Source  %
[noise_audio, Fs2] = audioread('AmbientNoise.wav');

%resample if different
[p,q] = rat(Fs/Fs2);
noise_audio = resample(noise_audio, p, q);
noise_audio = noise_audio(1:n);

%% positions of noise source, roughly parabolic with noise 2 feet behind and
% to sides 2x^2 - 2, mics at 0,0
noise_x = 2; % length is noise_x*2 + 1
noise_pos = horzcat((-noise_x:noise_x)',((-noise_x:noise_x)'.^2).*2 -2);


%% noise and audio overlaid
if speak
    sound(target_audio + noise_audio, Fs);
end

%% Step 2: Delay and attenuation due to distance from audio source
% attenuation by 1/(4pi*r)

%distance of each from target audio
mic_d = vecnorm((mic_pos - target_pos)');

% Subplot 1: Time-domain representation
figure;

subplot(mic_n+1, 1, 1);
plot(t, target_audio, 'b');

title('Time-Domain Representation of Audio Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

%delay of each and graph
SOUND = 343; %meters/second
mic_delay = mic_d/SOUND; 
mic_data = zeros(n, mic_n);
for i = 1:mic_n
    delay = int32(mic_delay(i)*Fs);

    %insert delay # of samples before and attenuate by 4pi*distance
    mic_data(:,i) = [zeros(delay,1); target_audio(1:end-delay)]/(4*pi*mic_d(i)) + add_noise(mic_pos(i), noise_pos, noise_audio, Fs);

    %adds noise to microphone
    %mic_data(:,i) = mic_data(:,i) + get_to_mic(mic_pos(i), noise_pos, noise_audio, Fs); 

    subplot(mic_n+1, 1, i+1);
    plot(t, mic_data(:,i));
    title(sprintf('Mic_{%d}, Delay %f ms', i, round(mic_delay(i)*1e3,3)));
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on; hold on;
end

%small test
%sound(mic_data(:,1),Fs);

if verbose
    fprintf("The microphones have delays of %s microseconds relative to mic1\n",join(string(abs(diff(mic_delay(1:2)))*1e6), ', ')); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       SIGNAL REACHES MICROPHONES                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% MICROPHONE CONTINUOUS TIME H(s)  

f_approx = [35, 100, 900, 1000, 1100, 8000, 15000]; %Hz
magnitude_dB = [-1.0, -0.5, 0.0, 0.0, 0.0, 1.5, 4.5];
phase_deg = [0, 0, 0, 0, 0, 0, 0,];
% Mic Frequency Response
[magnitudeCoeffs, phaseCoeffs] = fitFrequencyResponse(f_approx, magnitude_dB, phase_deg, 5);
mic_H = (10.^(evalFunction(magnitudeCoeffs,f)/20).*exp(1j*pi*evalFunction(phaseCoeffs,f)/180))';

%% ANTIALIASING FILTER, 6th order bessel filter to preserve phase
N=4;
fp = 2*pi*8000; %cutoff frequency
[z, p, k] = besself(N,fp);
[num, den] = zp2tf(z, p, k);
AAF = tf(num, den);
[H, ~] = freqresp(AAF,2*pi*f);
H = squeeze(H);

%%%% Plotting Group Delay and Magnitude Response of Filter
figure;
[h,w] = freqs(num,den);
subplot(2,1,1);
semilogx(w(1:end)/(2*pi),20*log(abs(h)));
title("Magnitude Response of AntiAliasing Filter")
xlabel("Frequency (Hz)")
ylabel("Magnitude (dB)")
xline(fp)
grid on

subplot(2,1,2);
grpdel = -diff(unwrap(angle(h)))./diff(w);
semilogx(w(2:end)/(2*pi),grpdel)
title("Group Delay of AntiAliasing Filter")
xlabel("Frequency (Hz)")
ylabel("Group delay (s)")
xline(fp)
grid on


%% %%%%%%%%% Apply Microphone and AAF tfs %%%%%%%%%%%%%%%%%%
input_H = H .* mic_H;
beam_input = fft(mic_data) .* repmat(input_H, 1,mic_n);

%%%%%%%%%%% FFT Graphs  %%%%%%%%%%%%%%%% Should Look V Similar
audio_fft = abs(fft(target_audio)); % Compute the FFT
beam_fft = abs(beam_input(:,1));

figure;
subplot(2,1,1);
semilogx(f(1:n/2), 20*log10(audio_fft(1:n/2,1))); % Plot only the positive frequencies
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
ylim([-200, 100]);
title('Magnitude Spectrum of Audio Signal');
subplot(2,1,2);

semilogx(f(1:n/2), 20*log10(beam_fft(1:n/2)), 'DisplayName','Beamforming Input'); hold on;
semilogx(f(1:n/2), 20*log10(abs(H(1:n/2))), 'b', 'LineWidth', 2,'DisplayName', 'AAF Filter');
semilogx(f(1:n/2), 20*log10(abs(mic_H(1:n/2))), 'r' , 'LineWidth', 2, 'DisplayName', 'Microphone H(s)');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
legend('Location', 'best');
ylim([-200, 100]);
title('Magnitude Spectrum of Audio Signal After Mics and AAF');

%%%%%%%%%%%%%%% FFT Check %%%%%%%%%%%%%

% SAMPLING - uses its own AAF, but doesnt matter if fs > 16kHz
fs_adc = 20000;
ts = (0:n-1)/fs_adc;
[p, q] = rat(fs_adc / Fs); % Calculate resampling factors

beam_input = ifft(beam_input);
beam_input = resample(beam_input, p, q); % Resample signal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         SIGNAL CONVERTED TO DIGITAL                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEAMFORMING OPERATIONS

% delay and sum --> simplest, least computationally intensive
  % After computing mic_data and mic_delay:
y_beamformed = delaysum(mic_data, mic_delay, Fs);
y_beamformed = y_beamformed / max(abs(y_beamformed));
audiowrite('beamformed_signal.wav', y_beamformed, Fs);

% Plot the beamformed signal in the time domain
figure;
subplot(2,1,1);
plot(t, y_beamformed, 'LineWidth', 1.5);
title('Beamformed Signal in Time Domain');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Compute the FFT of the beamformed signal for frequency-domain analysis
Y_Beamformed = fft(y_beamformed);
Y_Beamformed_magnitude = abs(Y_Beamformed);

% Plot the beamformed signal in the frequency domain (magnitude spectrum)
subplot(2,1,2);
semilogx(f(1:n/2), 20*log10(Y_Beamformed_magnitude(1:n/2)), 'LineWidth', 1.5);
title('Beamformed Signal Frequency Response');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;
ylim([-200, 100]);

% Plot time-domain: Overlay beamformed signal with the original
figure;
subplot(2,1,1);
plot(t, target_audio, 'DisplayName', 'Original', 'LineWidth', 1.5);
hold on;
plot(t, y_beamformed, 'DisplayName', 'Beamformed', 'LineWidth', 1.5);
hold off;
title('Time-Domain: Original vs Beamformed');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Location', 'best');
grid on;

% Frequency-domain comparison
Y_original = fft(target_audio);
Y_beamformed = fft(y_beamformed);

% Convert to magnitude (dB) for better comparison
Y_original_mag_dB = 20*log10(abs(Y_original(1:n/2)));
Y_beamformed_mag_dB = 20*log10(abs(Y_beamformed(1:n/2)));

subplot(2,1,2);
semilogx(f(1:n/2), Y_original_mag_dB, 'DisplayName', 'Original', 'LineWidth', 1.5);
hold on;
semilogx(f(1:n/2), Y_beamformed_mag_dB, 'DisplayName', 'Beamformed', 'LineWidth', 1.5);
hold off;
title('Frequency-Domain: Original vs Beamformed');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
legend('Location', 'best');
grid on;
ylim([-200, 100]);



% Bartlett Beamforming --> delay/sum + uniform weighting per microphone


% frequency domain delay/sum + deconvolution


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         COMPARISONS AND ANALYSIS                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% COMPARISON BETWEEN EACH


% fft audio data, multiply amplitude add phase

% Frequency points (Hz)




%%%%%%%                 HELPER FUNCTIONS                   %%%%%%%%%%%
function out = evalFunction(coeff, f)
    %from an array of coefficients evaluate function with log input
    out = 0;
    for i = 1:length(coeff)
        out = out + coeff(i)*log10(f).^(length(coeff)-i);
    end
end

% function takes source signal, sampling Fs, source locations, and mic position
% outputs resulting audio
function audio = add_noise(input_pos, sources, signal, Fs)
    d = vecnorm((sources - input_pos)'); %distance from mic to noise sources
    delays = int32((d/343)*Fs); 

    audio = zeros(length(signal),1);
    for i = 1:length(sources)
        audio = audio + [zeros(delays(i),1); signal(1:end-delays(i))]/(4*pi*d(i));
    end
    audio = audio/length(sources) .*2; 
end

