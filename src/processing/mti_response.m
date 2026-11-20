function H = mti_response(fd, PRI, order)
%MTI_RESPONSE  Amplitude response of an n-th order delay-line canceller (Lab 7, eqs. 7.2-7.3).
%   H = mti_response(fd, PRI, order)
%     fd     Doppler frequency [Hz] (any array)
%     PRI    pulse repetition interval [s]
%   |H_SC| = 2 |sin(omega PRI / 2)|,  |H_DC| = 4 sin^2(omega PRI / 2),  omega = 2 pi fd
%   i.e. |H| = (2 |sin(pi fd PRI)|)^order. The response is periodic in the PRF
%   and is zero at fd = k PRF: targets with those Dopplers are cancelled like
%   clutter (blind speeds v = k lambda PRF / 2).
H = (2*abs(sin(pi*fd*PRI))).^order;
end
