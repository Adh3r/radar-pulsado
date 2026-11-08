function k = radar_constants()
%RADAR_CONSTANTS  Physical constants used throughout the library.
%   k.c  speed of propagation [m/s] (academic value 3e8, as in the labs)
%   k.kB Boltzmann constant [J/K]
%   k.T0 reference noise temperature [K]
%   k.Z0 receiver input impedance [ohm]
k.c  = 3e8;
k.kB = 1.38e-23;
k.T0 = 290;
k.Z0 = 50;
end
