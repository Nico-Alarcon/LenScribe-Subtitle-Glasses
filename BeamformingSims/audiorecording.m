% Step 1: Set up the audio recorder
fs = 44100;          % Sample rate (Hz)
nBits = 16;          % Bit depth
nChannels = 1;       % Number of channels (1 = mono, 2 = stereo)
recObj = audiorecorder(fs, nBits, nChannels);
 
% Step 2: Record the audio
disp('Recording... Press any key to stop.');
record(recObj);      % Start recording
pause;               % Wait for user input to stop recording
stop(recObj);        % Stop recording
disp('Recording stopped.');

% Step 3: Retrieve and save the audio data
audioData = getaudiodata(recObj);               % Get the audio data
filename = 'recorded_audio.wav';                % Specify the filename
audiowrite(filename, audioData, fs);            % Save the audio to a .wav file

disp(['Audio saved as ', filename]);
