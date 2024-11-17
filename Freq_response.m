% Frequency points (Hz)
f = [20 40 60 80 100 200 400 600 800 1000 2000 4000 6000 8000 10000 20000];
% Amplitude response in dB (relative to a reference, e.g., 1kHz)
magnitude_dB = [-3 -1.25 -0.8 -0.40 -0.39 -0.05 0 0 0 0 0 0.01 0.18 0.2 1 4];
% Phase response (degrees)
phase_deg = [62 37 22 18 16 8 5 3 1 0.5 0 -1 -2 -4.8 -5.1 -5.2];
% Polynomial order
poly_order = 5;

% Call the function
[amplitudeFunc, phaseFunc] = fitFrequencyResponse(f, magnitude_dB, phase_deg, poly_order);

% Display the fitted functions
disp('Amplitude Response Function:');
disp(amplitudeFunc);

disp('Phase Response Function:');
disp(phaseFunc);
