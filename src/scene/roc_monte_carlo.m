function r = roc_monte_carlo(snr_chip_dB, code, n_trials, seed)
%ROC_MONTE_CARLO  Monte Carlo detection statistics with and without pulse compression.  [2026]
%   r = roc_monte_carlo(snr_chip_dB, code, n_trials, seed)
%   Each trial is a short fast-time record of unit-power complex noise plus a
%   non-fluctuating echo of the coded pulse (per-sample SNR = snr_chip_dB,
%   random phase). Square-law statistics are collected
%     * 'raw'        : |x|^2 at one sample of the echo, no compression
%     * 'compressed' : |y|^2 at the matched-filter peak, normalised so that
%                      the output noise power is 1
%   for H1 (target cell) and H0 (cells away from the target).
%   r.h1_raw, r.h0_raw, r.h1_mf, r.h0_mf  detector statistics
%   r.gain_dB  theoretical compression gain 10 log10(numel(code))
rng(seed);
code = code(:);
K  = numel(code);
nf = 8*K;                         % record length
n0 = 3*K;                         % echo start
cells0 = [1:(n0-K), (n0+2*K):nf]; % H0 cells, clear of the echo and its sidelobes

amp = sqrt(db2lin(snr_chip_dB));
x = sqrt(1/2)*(randn(nf, n_trials) + 1i*randn(nf, n_trials));
x(n0:n0+K-1, :) = x(n0:n0+K-1, :) + amp * code .* exp(1i*2*pi*rand(1, n_trials));

y = matched_filter(x, code, 'aligned') / sqrt(K);   % unit output noise power

r.h1_raw = abs(x(n0, :)).^2;
r.h0_raw = abs(x(cells0, :)).^2;
r.h1_mf  = abs(y(n0, :)).^2;
r.h0_mf  = abs(y(cells0, :)).^2;
r.gain_dB = 10*log10(K);
end
