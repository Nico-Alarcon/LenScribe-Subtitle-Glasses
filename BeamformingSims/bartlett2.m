function [power_output, y] = bartlett2(x,d,N,theta)
    Rxx = (x' * x)/size(x,2);
    steering_vector = exp(-1j * 2 * pi * d * (0:N-1) * sin(theta));
    power_output = steering_vector * Rxx * steering_vector';
    y = x * steering_vector'/N;
end