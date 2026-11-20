function sigma = rcs_circular_plate(r, freq, theta_deg)
%RCS_CIRCULAR_PLATE  Physical-optics monostatic RCS of a PEC circular flat plate [m^2].
%   sigma = rcs_circular_plate(r, freq, theta_deg)
%     r          plate radius [m]
%     freq       frequency [Hz]
%     theta_deg  aspect angle from the plate normal [deg] (any array)
%
%     sigma(theta) = (4 pi^3 r^4 / lambda^2) * [ 2 J1(x)/x * cos(theta) ]^2,
%     x = 2 k r sin(theta) = (4 pi r / lambda) sin(theta)
%
%   At normal incidence it reduces to 4 pi A^2 / lambda^2 with A = pi r^2.
%   Physical optics ignores edge diffraction, so it goes to zero edge-on
%   (theta = 90 deg), where the real RCS is small but not zero.
k = radar_constants();
lambda = k.c/freq;
th = deg2rad(theta_deg);
x  = 4*pi*r/lambda * sin(th);
airy = ones(size(x));                       % 2 J1(x)/x -> 1 as x -> 0
nz = abs(x) > 1e-12;
airy(nz) = 2*besselj(1, x(nz)) ./ x(nz);
sigma = 4*pi^3*r^4/lambda^2 * (airy .* cos(th)).^2;
end
