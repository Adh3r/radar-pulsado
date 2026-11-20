function [A, delay, fd] = ambiguity_function(s, Ts, fd)
%AMBIGUITY_FUNCTION  Normalised range-Doppler ambiguity |chi(tau, fd)| of a pulse.  [2026]
%   [A, delay, fd] = ambiguity_function(s, Ts, fd)
%     s   baseband pulse samples, Ts sampling time [s], fd Doppler axis [Hz]
%   A(i, :) is the matched-filter output (filter matched to the zero-Doppler
%   pulse) for an echo shifted by fd(i), normalised to the zero-Doppler peak.
%   delay is the lag axis [s]. Row fd = 0 is the autocorrelation; the
%   fd-axis cut shows how fast a Doppler mismatch erodes the peak.
s = s(:).';
N = numel(s);
A = zeros(numel(fd), 2*N - 1);
for i = 1:numel(fd)
    A(i, :) = abs(matched_filter(apply_doppler(s, fd(i), Ts), s));
end
A = A / sum(abs(s).^2);
delay = (-(N-1):(N-1)) * Ts;
end
