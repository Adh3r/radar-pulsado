function [m, v, t, env] = pulsed_radar(tau, PRI, fc, Pt_dBW, L_dB, Np, fs, R)
%PULSED_RADAR  Lab 1.1 interface: radar figures plus the received pulse burst.
%   [m, v, t, env] = pulsed_radar(tau, PRI, fc, Pt_dBW, L_dB, Np, fs, R)
%   m is the struct returned by PULSE_METRICS; v, t, env as in PULSE_TRAIN.
m = pulse_metrics(tau, PRI, Pt_dBW, L_dB);
[v, t, env] = pulse_train(tau, PRI, fc, Pt_dBW, L_dB, Np, fs, R);
end
