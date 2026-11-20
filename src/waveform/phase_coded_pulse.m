function [s, noise, sn, Ts] = phase_coded_pulse(tau_pulse, n_samp_chip, phases_deg, amps, SNR_dB)
%PHASE_CODED_PULSE  Baseband phase-coded pulse plus complex white Gaussian noise (Lab 4.1).
%   [s, noise, sn, Ts] = phase_coded_pulse(tau_pulse, n_samp_chip, phases_deg, amps, SNR_dB)
%     tau_pulse   total pulse duration [s]
%     n_samp_chip samples per chip
%     phases_deg  chip phases [deg] (Barker: 0 -> +1, 180 -> -1)
%     amps        chip amplitudes ([] -> all ones)
%     SNR_dB      per-sample SNR = mean |s|^2 / noise power
%   Outputs: clean signal s, noise, noisy signal sn = s + noise, sampling time Ts.
%   Noise is drawn with randn: call rng(...) beforehand for reproducible runs.
K = numel(phases_deg);
if nargin < 4 || isempty(amps), amps = ones(1, K); end

Ts = tau_pulse / K / n_samp_chip;
s  = repelem(amps(:).' .* exp(1i*deg2rad(phases_deg(:).')), n_samp_chip);

Pn    = mean(abs(s).^2) / db2lin(SNR_dB);
noise = sqrt(Pn/2) * (randn(size(s)) + 1i*randn(size(s)));
sn    = s + noise;
end
