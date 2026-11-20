function c = clutter_slow_time(n_cells, n_pulses, sigma_f, PRI, power)
%CLUTTER_SLOW_TIME  Complex Gaussian clutter with a Gaussian Doppler spectrum (eq. 7.1).  [2026]
%   c = clutter_slow_time(n_cells, n_pulses, sigma_f, PRI, power)
%     n_cells   independent range cells (rows)
%     n_pulses  pulses (columns)
%     sigma_f   spectral spread [Hz] (0 -> perfectly static clutter)
%     power     clutter power per cell: scalar or n_cells x 1
%   Each row is a zero-mean complex Gaussian sequence with autocorrelation
%   rho(m) = exp(-2 pi^2 sigma_f^2 (m PRI)^2), centred at zero Doppler. It is
%   built by colouring white noise with a square root of the Toeplitz
%   correlation matrix (eigen-decomposition, robust when rho is close to 1).
m = 0:n_pulses-1;
R = toeplitz(exp(-2*pi^2*sigma_f^2*(m*PRI).^2));
[V, D] = eig((R + R')/2);
A = V * diag(sqrt(max(diag(D), 0)));          % R = A A'
w = sqrt(1/2)*(randn(n_cells, n_pulses) + 1i*randn(n_cells, n_pulses));
c = sqrt(power(:)) .* (w * A.');
end
