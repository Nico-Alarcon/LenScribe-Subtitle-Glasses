% Frequency points (Hz)
f = [35, 100, 900, 1000, 1100, 8000, 15000];
% Amplitude response in dB (relative to a reference, e.g., 1kHz)
magnitude_dB = [-1.0, -0.5, 0.0, 0.0, 0.0, 1.5, 4.5];
% Phase response (degrees)
phase_deg = [0, 0, 0, 0, 0, 0, 0,];
% Polynomial order
poly_order = 5;

% Call the function
[magnitudeCoeffs, phaseCoeffs] = fitFrequencyResponse(f, magnitude_dB, phase_deg, poly_order);

% Display the coefficients
disp('Magnitude Response Coefficients:');
disp(magnitudeCoeffs);

disp('Phase Response Coefficients:');
disp(phaseCoeffs);
