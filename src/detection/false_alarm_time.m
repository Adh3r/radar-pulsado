function [Tfa, Pfa] = false_alarm_time(ratio_dB, BW)
%FALSE_ALARM_TIME  Mean time between false alarms of an envelope detector (Lab 3.1).
%   [Tfa, Pfa] = false_alarm_time(ratio_dB, BW)
%     ratio_dB  threshold-to-noise ratio VT^2/(2 psi0) [dB] (vector)
%     BW        receiver bandwidth [Hz] (vector)
%   Tfa(i,j) [s] for ratio_dB(i) and BW(j); Pfa(i) per threshold.
%
%   Gaussian noise -> Rayleigh envelope: Pfa = exp(-VT^2/(2 psi0)).
%   The receiver offers about BW independent decisions per second, so
%       Tfa = 1/(BW Pfa) = exp(VT^2/(2 psi0)) / BW
%   Exponential in the threshold, inverse in the bandwidth.
ratio_lin = db2lin(ratio_dB(:));       % power ratio -> /10
Pfa = exp(-ratio_lin);
Tfa = (1./Pfa) ./ BW(:).';
end
