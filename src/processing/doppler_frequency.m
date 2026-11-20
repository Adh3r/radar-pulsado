function fd = doppler_frequency(vr, f0)
%DOPPLER_FREQUENCY  Two-way Doppler shift fd = 2 vr f0 / c [Hz] (vr > 0: approaching).
%   The Lab 4 statement writes omega_d = 2 f0 2 vr / c without the 2*pi;
%   the physically correct angular frequency is 2*pi*fd.
k = radar_constants();
fd = 2*vr.*f0/k.c;
end
