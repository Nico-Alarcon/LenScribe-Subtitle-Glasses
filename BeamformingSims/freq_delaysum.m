function y = freq_delaysum(X, f, delays)
    steering_matrix = exp(-1j * 2 * pi * f' * delays);
    y = mean(X .* steering_matrix,2);
end