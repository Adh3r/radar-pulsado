function sd = apply_doppler(s, fd, Ts)
%APPLY_DOPPLER  Impose a Doppler shift on baseband samples: s(t) exp(j 2 pi fd t).
%   sd = apply_doppler(s, fd, Ts), t = (0:N-1) Ts along the vector s.
t  = (0:numel(s)-1) * Ts;
sd = s .* reshape(exp(1i*2*pi*fd*t), size(s));
end
