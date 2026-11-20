function [CA_dB, I_dB] = mti_clutter_attenuation(sigma_f, PRI, order)
%MTI_CLUTTER_ATTENUATION  Clutter attenuation and improvement factor of a canceller.  [2026]
%   [CA_dB, I_dB] = mti_clutter_attenuation(sigma_f, PRI, order)
%     sigma_f  standard deviation of the Gaussian clutter spectrum [Hz]
%              (eq. 7.1 with omega_0 = 0; sigma_omega = 2 pi sigma_f)
%   A Gaussian spectrum has autocorrelation rho(tau) = exp(-2 pi^2 sigma_f^2 tau^2),
%   so for a canceller with taps h the output clutter power is
%       Pc_out / Pc_in = sum_i sum_j h_i h_j rho((i - j) PRI)
%   CA = Pc_in / Pc_out (the lab's figure of merit). The improvement factor
%   I = CA * sum(h^2) also credits the average gain on targets whose Doppler
%   is uniformly distributed (2 for the single, 6 for the double canceller).
k = 0:order;
h = (-1).^k .* arrayfun(@(i) nchoosek(order, i), k);
[I, J] = ndgrid(k, k);
lag = abs(I - J) * PRI;
CA_dB = zeros(size(sigma_f));
for n = 1:numel(sigma_f)
    rho = exp(-2*pi^2*sigma_f(n)^2*lag.^2);
    CA_dB(n) = -10*log10(h * rho * h.');
end
I_dB = CA_dB + 10*log10(sum(h.^2));
end
