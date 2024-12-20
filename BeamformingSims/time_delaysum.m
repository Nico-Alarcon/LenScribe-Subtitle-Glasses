function y = time_delaysum(X, delays, fs)

    % Input validation
    [numSamples, numMics] = size(X);
    if length(delays) ~= numMics
        error('Number of delays must match number of microphone channels in X.');
    end
    
    % Convert delays to sample indices
    sampleShifts = delays * fs;
    t = (0:numSamples-1)'; 
    y_delayed = zeros(numSamples, numMics);
    
    for m = 1:numMics
        shift = sampleShifts(m);
        
        if abs(shift) < 1e-12 %no delay needed
            y_delayed(:, m) = X(:, m);
        else
            t_shifted = t - shift; % Fractional delay via linear interpolation

            t_shifted_clamped = min(max(t_shifted, 0), numSamples - 1);

            idx_floor = floor(t_shifted_clamped);
            frac = t_shifted_clamped - idx_floor;
            idx_ceil = idx_floor + 1;
            idx_ceil(idx_ceil > (numSamples-1)) = numSamples - 1;

            x_floor = X(idx_floor+1, m); % Interpolate 
            x_ceil = X(idx_ceil+1, m);
            y_delayed(:, m) = (1 - frac).*x_floor + frac.*x_ceil;
        end
    end

    % Average across microphones
    y = mean(y_delayed, 2);
end
