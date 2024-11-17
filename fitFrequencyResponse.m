function [amplitudeFunc, phaseFunc] = fitFrequencyResponse(f, magnitude_dB, phase_deg, poly_order)
    % FITFREQUENCYRESPONSE Fits amplitude and phase responses to polynomial functions.
    % 
    % INPUTS:
    % f - Frequency points (Hz)
    % magnitude_dB - Amplitude response in dB
    % phase_deg - Phase response in degrees
    % poly_order - Order of polynomial for fitting
    %
    % OUTPUTS:
    % amplitudeFunc - String representation of fitted amplitude function
    % phaseFunc - String representation of fitted phase function

    % Define finer frequency points for interpolation
    f_interp = logspace(log10(min(f)), log10(max(f)), 1000);

    % Interpolate the amplitude response
    magnitude_interp = interp1(f, magnitude_dB, f_interp, 'pchip');

    % Interpolate the phase response
    phase_interp = interp1(f, phase_deg, f_interp, 'pchip');

    % Fit polynomial functions
    magnitude_fit_coeffs = polyfit(log10(f_interp), magnitude_interp, poly_order);
    phase_fit_coeffs = polyfit(log10(f_interp), phase_interp, poly_order);

    % Construct the amplitude function as a string
    amplitudeFunc = 'y = ';
    for i = 1:length(magnitude_fit_coeffs)
        term = sprintf('%.4f * log10(f)^%d', magnitude_fit_coeffs(i), poly_order - i + 1);
        amplitudeFunc = strcat(amplitudeFunc, term);
        if i < length(magnitude_fit_coeffs)
            amplitudeFunc = strcat(amplitudeFunc, ' + ');
        end
    end

    % Construct the phase function as a string
    phaseFunc = 'y = ';
    for i = 1:length(phase_fit_coeffs)
        term = sprintf('%.4f * log10(f)^%d', phase_fit_coeffs(i), poly_order - i + 1);
        phaseFunc = strcat(phaseFunc, term);
        if i < length(phase_fit_coeffs)
            phaseFunc = strcat(phaseFunc, ' + ');
        end
    end
end
