function T = threshold_from_pfa(Pfa, noise_power, n_pulses)
%THRESHOLD_FROM_PFA  Square-law detection threshold for a target false-alarm probability.  [2026]
%   T = threshold_from_pfa(Pfa, noise_power, n_pulses)
%     Pfa          desired false-alarm probability per cell
%     noise_power  mean noise power sigma^2 at the detector input (|y|^2 domain)
%     n_pulses     pulses integrated non-coherently (sum of |y|^2), default 1
%   With complex Gaussian noise each |y|^2 is exponential with mean sigma^2;
%   the sum of n of them is Gamma(n, sigma^2), hence
%       T = sigma^2 * gammaincinv(Pfa, n, 'upper')   (n = 1: T = -sigma^2 ln Pfa)
%   This is a FIXED threshold: it assumes the noise power is known. CFAR
%   (CA_CFAR, Lab 6) estimates it locally instead.
if nargin < 3, n_pulses = 1; end
T = noise_power .* gammaincinv(Pfa, n_pulses, 'upper');
end
