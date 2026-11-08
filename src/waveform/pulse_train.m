function [v, t, env] = pulse_train(tau, PRI, fc, Pt_dBW, L_dB, Np, fs, R)
%PULSE_TRAIN  Received RF pulse burst from a target at range R (steady state).
%   [v, t, env] = pulse_train(tau, PRI, fc, Pt_dBW, L_dB, Np, fs, R)
%     Np pulses observed, sampled at fs. v is the RF voltage on Z0 [V], t the
%     time axis [s] and env the 0/1 envelope.
%   The radar is assumed to have been transmitting before t = 0, so echoes
%   of pulses sent before the window appear inside it: when 2R/c > PRI the
%   echo lands at the folded (ambiguous) delay, as a real receiver shows it.
k  = radar_constants();
m  = pulse_metrics(tau, PRI, Pt_dBW, L_dB);
N  = round(Np*PRI*fs);
t  = (0:N-1)/fs;
td = 2*R/k.c;                                    % round-trip delay

% Pulses are placed on integer sample indices: every pulse has exactly
% M = round(tau*fs) samples. (Comparing floating-point times, t >= t0 &
% t < t0 + tau, gave 100 or 101 samples depending on rounding: a bug in the
% original lab code, found by the closed-form spectrum check in s01.)
M   = round(tau*fs);
env = zeros(1, N);
for p = -ceil(td/PRI):(Np-1)                     % include pulses sent before t = 0
    n0 = ceil((p*PRI + td)*fs - 1e-6);           % first sample at or after the echo start
    idx = n0 + (1:M);                            % MATLAB indices (sample n -> index n+1)
    env(idx(idx >= 1 & idx <= N)) = 1;
end
v = m.Vpeak * env .* cos(2*pi*fc*t);
end
