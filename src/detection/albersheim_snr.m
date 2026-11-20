function [SNR_dB, SNR_lin] = albersheim_snr(Pd, Pfa)
%ALBERSHEIM_SNR  Single-pulse SNR required for (Pd, Pfa), Albersheim's equation (Lab 3.5).
%   [SNR_dB, SNR_lin] = albersheim_snr(Pd, Pfa)
%   SNR_dB(i,j) for Pd(i) and Pfa(j). Empirical fit, valid for
%   0.3 <= Pd <= 0.99; its error against the exact Marcum-Q result is
%   measured in scripts/s03_detection_theory.m and tests/test_detection.m.
%
%       SNR = A + 0.12 A B + 1.7 B,   A = ln(0.62/Pfa),  B = ln(Pd/(1-Pd))
%
%   Note: the course formula sheet printed the last coefficient as 0.17;
%   1.7 is the value that reproduces Skolnik's curves (Pd = 0.9, Pfa = 1e-6
%   -> 13.1 dB).
if any(Pd(:) < 0.3 | Pd(:) > 0.99)
    warning('albersheim_snr:range', 'Albersheim is only valid for 0.3 <= Pd <= 0.99.');
end
A = log(0.62 ./ Pfa(:).');        % row: stricter false-alarm -> larger A
B = log(Pd(:) ./ (1 - Pd(:)));    % column: logit of Pd, explodes as Pd -> 1
SNR_lin = A + 0.12*B.*A + 1.7*B;
SNR_dB  = lin2db(SNR_lin);
end
