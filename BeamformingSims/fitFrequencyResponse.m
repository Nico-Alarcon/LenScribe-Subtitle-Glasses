function [magnitudeCoeffs, phaseCoeffs] = fitFrequencyResponse(f, magnitude_dB, phase_deg, poly_order)
    % FITFREQUENCYRESPONSE Fits amplitude and phase responses to polynomial functions.
    % 
    % INPUTS:
    % f - Frequency points (Hz)
    % magnitude_dB - Amplitude response in dB
    % phase_deg - Phase response in degrees
    % poly_order - Order of polynomial for fitting
    %
    % OUTPUTS:
    % magnitudeCoeffs - Coefficients of the fitted polynomial for magnitude response
    % phaseCoeffs - Coefficients of the fitted polynomial for phase response

    % Define finer frequency points for interpolation
    f_interp = logspace(log10(min(f)), log10(max(f)), 1000);

    % Interpolate the amplitude response
    magnitude_interp = interp1(f, magnitude_dB, f_interp, 'pchip');

    % Interpolate the phase response
    phase_interp = interp1(f, phase_deg, f_interp, 'pchip');

    % Fit polynomial functions
    magnitudeCoeffs = polyfit(log10(f_interp), magnitude_interp, poly_order);
    phaseCoeffs = polyfit(log10(f_interp), phase_interp, poly_order);
end
