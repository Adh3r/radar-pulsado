function sigma_dBsm = rcs_scatterers(sigma_n_dBsm, pos_xyz, radar_sph, fc)
%RCS_SCATTERERS  Monostatic RCS of an object modelled as a set of point scatterers (Lab 5.3).
%   sigma_dBsm = rcs_scatterers(sigma_n_dBsm, pos_xyz, radar_sph, fc)
%     sigma_n_dBsm  column vector of scatterer (sphere) RCSs [dBsm]
%     pos_xyz       N x 3 scatterer positions [x y z] [m]
%     radar_sph     M x 3 radar positions in spherical coordinates
%                   [r [m], theta [deg], phi [deg]]: theta is the polar angle
%                   from +z (theta = 90 deg is the XY plane), phi the azimuth from +x
%     fc            carrier frequency [Hz]
%   sigma_dBsm      M x 1 monostatic RCS seen from each radar position [dBsm]
%
%   The field of each scatterer is sqrt(sigma_n) with the two-way phase
%   2 k r_n to that radar position (exact distances, no far-field
%   approximation), and the contributions add coherently:
%       sigma = | sum_n sqrt(sigma_n) exp(-j 2 k r_n) |^2
k = radar_constants();
kw = 2*pi*fc/k.c;
r  = radar_sph(:, 1); th = deg2rad(radar_sph(:, 2)); ph = deg2rad(radar_sph(:, 3));
radar_xyz = [r.*sin(th).*cos(ph), r.*sin(th).*sin(ph), r.*cos(th)];          % M x 3
amp = sqrt(db2lin(sigma_n_dBsm(:))).';                                        % 1 x N
d = sqrt((radar_xyz(:,1) - pos_xyz(:,1).').^2 + ...
         (radar_xyz(:,2) - pos_xyz(:,2).').^2 + ...
         (radar_xyz(:,3) - pos_xyz(:,3).').^2);                               % M x N distances
E = sum(amp .* exp(-1i*2*kw*d), 2);
sigma_dBsm = lin2db(abs(E).^2);
end
