function m = pulse_metrics(tau, PRI, Pt_dBW, L_dB)
%PULSE_METRICS  Range/timing and power figures of a simple (uncoded) pulsed radar.
%   m = pulse_metrics(tau, PRI, Pt_dBW, L_dB)
%     tau    pulse width [s]           PRI   pulse repetition interval [s]
%     Pt_dBW peak transmitted power    L_dB  link losses [dB]
%   Fields of m:
%     PRF [Hz], Rmin blind range [m], Rmax unambiguous range [m], dR range
%     resolution [m], duty [%], Ppeak_rx [W], Pavg_dBm received mean power,
%     Epulse received pulse energy [J], Vpeak peak voltage on Z0 [V]
%
%   Every range figure follows from the round-trip relation R = c*t/2.
k = radar_constants();
m.PRF  = 1/PRI;
m.Rmin = k.c*tau/2;          % receiver blanked while transmitting
m.Rmax = k.c*PRI/2;          % echo must return before the next pulse
m.dR   = k.c*tau/2;          % uncompressed pulse: resolution = c*tau/2
m.duty = 100*tau/PRI;

m.Ppeak_rx = db2lin(Pt_dBW - L_dB);              % [W]
m.Pavg_dBm = lin2db(m.Ppeak_rx*tau/PRI*1e3);     % mean = peak * duty cycle
m.Epulse   = m.Ppeak_rx*tau;                     % [J]
m.Vpeak    = sqrt(2*m.Ppeak_rx*k.Z0);            % sinusoid: P = Vpeak^2/(2*Z0)
end
