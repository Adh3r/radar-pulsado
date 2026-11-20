%% S05 - RADAR CROSS SECTION: simple shapes, multi-scatterer targets, Swerling models   (Lab 5)
%  Exercises 5.1-5.6 using src/rcs. The sphere and plate are own
%  implementations of the exact Mie series and of physical optics (the
%  course handout provides adapted routines from Mahafza, 2005, which are not
%  redistributed here; they were used locally to cross-check these results).
%  Extension 2026: what RCS fluctuation costs in detection (Swerling 0/I/III).

clear; close all;
run(fullfile(fileparts(mfilename('fullpath')), '..', 'startup.m'));
rng(5);
k = radar_constants();

%% 5.1 Sphere: frequency dependence (Rayleigh, Mie, optical regions)
ka = 0.05:0.01:15;                               % sphere circumference in wavelengths
s = rcs_sphere_mie(ka);                          % sigma / (pi a^2)
[pk, ipk] = max(s);
fprintf('\n=== 5.1 PEC sphere ===\n');
fprintf('Rayleigh check, ka = 0.05: %.4e vs 9 (ka)^4 = %.4e\n', s(1), 9*ka(1)^4);
fprintf('Mie: first resonance peak %.2f at ka = %.2f\n', pk, ka(ipk));
fprintf('optical: mean sigma/(pi a^2) for ka in [10, 15] = %.3f (-> 1, i.e. sigma = pi a^2)\n', mean(s(ka >= 10)));

fig = figure('Name', '5.1 Sphere RCS', 'Position', [100 100 1100 400]);
subplot(1,2,1);
plot(ka, s, 'b', 'LineWidth', 1.3); grid on; hold on; yline(1, 'k--', '\pi a^2');
xlabel('Circumference in wavelengths, 2\pi a/\lambda'); ylabel('\sigma / \pi a^2');
title('Normalised RCS of a PEC sphere (linear)'); xlim([0 15]);
subplot(1,2,2);
loglog(ka, s, 'b', 'LineWidth', 1.3); grid on; hold on;
loglog(ka(ka < 1), 9*ka(ka < 1).^4, 'r--', 'LineWidth', 1.3);
loglog([1 15], [1 1], 'k--');
xline(1, 'k:'); xline(10, 'k:');
text(0.12, 2, 'Rayleigh'); text(2.2, 6, 'Mie (resonance)'); text(10.5, 2.5, 'optical');
xlabel('2\pi a/\lambda'); ylabel('\sigma / \pi a^2'); ylim([1e-5 10]);
title('Log scale: \sigma \propto f^4, then oscillation, then \pi a^2');
legend('Mie series', '9 (ka)^4 (Rayleigh)', '\pi a^2 (optical)', 'Location', 'southeast');
save_figure(fig, 's05_sphere_rcs');
% Discussion: for a sphere much smaller than lambda (Rayleigh), sigma grows
% as f^4 (40 dB/decade). Around ka ~ 1 the specular return and the creeping
% wave around the sphere interfere (Mie/resonance region, peak ~3.6 pi a^2).
% For ka >> 1 (optical region) the creeping wave dies out and sigma -> pi a^2,
% independent of frequency: the projected area.
% Note: the handout routine returns |sum| rather than |sum|^2, i.e. the square
% root of the normalised RCS, so its dB curve shows half the physical slope
% (20 instead of 40 dB/decade in the Rayleigh region).

%% 5.2 Circular flat plate: aspect-angle dependence
r = 0.25; f_p = 3e9; lambda_p = k.c/f_p;          % values chosen here (2.5 lambda radius)
th = 0:0.05:90;
sp = rcs_circular_plate(r, f_p, th);
A = pi*r^2;
fprintf('\n=== 5.2 Circular plate (r = %.2f m, %.0f GHz) ===\n', r, f_p/1e9);
fprintf('normal incidence: %.2f m^2 (%.2f dBsm); 4 pi A^2/lambda^2 = %.2f m^2\n', sp(1), lin2db(sp(1)), 4*pi*A^2/lambda_p^2);
fprintf('first null at %.2f deg (theory asin(3.8317 lambda / (4 pi r)) = %.2f deg)\n', ...
        th(find(diff(sign(diff(sp))) > 0, 1) + 1), asind(3.8317*lambda_p/(4*pi*r)));
fig = figure('Name', '5.2 Plate RCS');
plot(th, lin2db(sp + 1e-12), 'b', 'LineWidth', 1.3); grid on;
yline(lin2db(4*pi*A^2/lambda_p^2), 'r--', '4\piA^2/\lambda^2');
xlabel('Aspect angle from the normal [deg]'); ylabel('RCS [dBsm]'); ylim([-60 25]);
title(sprintf('PEC circular plate, r = %.2f m = %.1f\\lambda, %.0f GHz (physical optics)', r, r/lambda_p, f_p/1e9));
save_figure(fig, 's05_plate_rcs');
% Discussion: at normal incidence (0 deg) the whole plate is equidistant from
% the radar, all its surface currents re-radiate in phase and the plate acts
% like an antenna of aperture A and gain 4 pi A / lambda^2, hence
% 4 pi A^2 / lambda^2. Tilting it adds a linear phase across the surface, the
% contributions cancel (Airy pattern 2 J1(x)/x, first null at ~7 deg here),
% and edge-on (90 deg) the projected area is zero. Physical optics predicts
% exactly zero there; the real edge diffraction leaves a small residual RCS.

%% 5.4 Two 1 m^2 spheres at (-L/2, 0, 0) and (L/2, 0, 0), 1.5 GHz
fc4 = 1.5e9; lambda4 = k.c/fc4; kw = 2*pi/lambda4;
phi = (0:0.1:360).';
radar = [1e3*ones(size(phi)), 90*ones(size(phi)), phi];          % 1 km, XY plane
Ls = [1 0.1 0.01];
fig = figure('Name', '5.4 Two spheres', 'Position', [60 60 1200 700]);
fprintf('\n=== 5.4 Two 1 m^2 spheres, %.1f GHz (lambda = %.0f cm) ===\n', fc4/1e9, lambda4*100);
for i = 1:3
    L = Ls(i);
    sig = rcs_scatterers([0; 0], [-L/2 0 0; L/2 0 0], radar, fc4);
    sig_th = lin2db(4*cos(kw*L*cosd(phi)).^2 + 1e-12);           % far-field closed form
    lin = db2lin(sig);
    lc = lin(1:end-1);                                           % 0..359.9 deg, circular
    nmax = nnz(lc > circshift(lc, 1) & lc >= circshift(lc, -1) & lc > 3.9);
    fprintf('L = %5.2f m (%.2f lambda): max %.2f dBsm, %d maxima over 360 deg, max |sim - far field| = %.3f dB (where > -20 dBsm)\n', ...
            L, L/lambda4, max(sig), nmax, max(abs(sig - sig_th) .* (sig_th > -20)));
    subplot(2,3,i);
    plot(phi, sig, 'b', 'LineWidth', 1.1); grid on; xlim([0 360]); ylim([-30 10]);
    xlabel('Aspect angle \phi [deg]'); ylabel('RCS [dBsm]');
    title(sprintf('L = %g m = %.2f\\lambda', L, L/lambda4));
    ax = subplot(2,3,3+i);
    pax = polaraxes('Position', ax.Position); delete(ax);
    polarplot(pax, deg2rad(phi), lin, 'b', 'LineWidth', 1.1);
    pax.RLim = [0 4.2]; title(pax, 'linear [m^2]');
end
save_figure(fig, 's05_two_spheres');
% Discussion: the two echoes add with a path difference 2 L cos(phi), so
% sigma = 4 sigma0 cos^2(k L cos(phi)): between 0 (opposite phase) and
% 4 sigma0 = +6 dBsm (in phase, the coherent sum of two 1 m^2 spheres). The
% larger the electrical separation L/lambda, the faster the phase changes with
% aspect: more and narrower lobes (L = 1 m = 5 lambda). At L = 1 cm =
% lambda/20 the two spheres are never out of phase and behave as a single
% 4 m^2 scatterer, practically isotropic.

%% 5.5-5.6 Swerling statistics from 10 random scatterers, 3 GHz, 1 km
fc5 = 3e9; R_circle = 5;                          % "circle of 10 m": 10 m diameter assumed
phi5 = (0:0.01:360).';                            % steps of 0.01 deg
radar5 = [1e3*ones(size(phi5)), 90*ones(size(phi5)), phi5];
ns = 10;
rr = R_circle*sqrt(rand(ns,1)); aa = 2*pi*rand(ns,1);
pos = [rr.*cos(aa), rr.*sin(aa), zeros(ns,1)];   % uniform in the disc
sig_n = 0.5*randn(ns, 1);                         % RCSs close to 0 dBsm
sw12 = db2lin(rcs_scatterers(sig_n, pos, radar5, fc5));
sig_n3 = sig_n; sig_n3(1) = sig_n3(1) + 17;       % one scatterer 17 dB above the others
sw34 = db2lin(rcs_scatterers(sig_n3, pos, radar5, fc5));

u12 = sw12/mean(sw12); u34 = sw34/mean(sw34);
x = linspace(0, 6, 400);
Kr = db2lin(sig_n3(1)) / sum(db2lin(sig_n3(2:end)));         % dominant / rest power ratio
Pr = 1/(1 + Kr); Pd_ = Kr/(1 + Kr);                          % normalised to unit mean
rice = exp(-(x + Pd_)/Pr) .* besseli(0, 2*sqrt(x*Pd_)/Pr) / Pr;   % dominant + Gaussian rest
ks = @(u, cdf) max(abs((1:numel(u))'/numel(u) - cdf(sort(u))));
cdf_sw1 = @(v) 1 - exp(-v);
cdf_sw3 = @(v) 1 - (1 + 2*v).*exp(-2*v);
cdf_rice = @(v) interp1(x, cumtrapz(x, rice), min(v, x(end)));
fprintf('\n=== 5.5-5.6 Swerling statistics (10 scatterers, %d aspect samples) ===\n', numel(phi5));
fprintf('5.5 equal scatterers : std/mean = %.2f (exponential: 1); KS distance to SW I/II %.3f, to SW III/IV %.3f\n', ...
        std(u12), ks(u12, cdf_sw1), ks(u12, cdf_sw3));
fprintf('5.6 one +17 dB       : std/mean = %.2f (chi-square 4: 0.71); KS to SW I/II %.3f, to SW III/IV %.3f, to Rice (K = %.1f) %.3f\n', ...
        std(u34), ks(u34, cdf_sw1), ks(u34, cdf_sw3), Kr, ks(u34, cdf_rice));

fig = figure('Name', '5.5-5.6 Swerling', 'Position', [100 100 1100 420]);
subplot(1,2,1);
histogram(u12, 80, 'Normalization', 'pdf', 'FaceColor', [0.6 0.7 0.9]); hold on; grid on;
plot(x, swerling_pdf(x, 1, 1), 'r', 'LineWidth', 1.6); plot(x, swerling_pdf(x, 1, 3), 'k--', 'LineWidth', 1.3);
xlim([0 6]); xlabel('\sigma / mean(\sigma)'); ylabel('pdf');
title('5.5: 10 comparable scatterers'); legend('simulated', 'SW I/II (exponential)', 'SW III/IV (\chi^2_4)');
subplot(1,2,2);
histogram(u34, 80, 'Normalization', 'pdf', 'FaceColor', [0.6 0.7 0.9]); hold on; grid on;
plot(x, swerling_pdf(x, 1, 3), 'r', 'LineWidth', 1.6); plot(x, swerling_pdf(x, 1, 1), 'k--', 'LineWidth', 1.3);
plot(x, rice, 'm:', 'LineWidth', 1.8);
xlim([0 4]); xlabel('\sigma / mean(\sigma)');
title('5.6: one scatterer 17 dB above the rest');
legend('simulated', 'SW III/IV (\chi^2_4)', 'SW I/II (exponential)', sprintf('dominant + Rayleigh (K = %.1f)', Kr));
save_figure(fig, 's05_swerling_statistics');
% Discussion: with comparable scatterers the RCS is the squared magnitude of
% a sum of random phasors, a complex Gaussian: exponential pdf (std = mean),
% the Swerling I/II model. With one dominant scatterer the pdf loses its
% peak at zero (deep fades become rare): that is what Swerling III/IV
% (chi-square with 4 dof) models. The chi-square-4 law is itself an
% approximation of "dominant + Rayleigh"; with a dominant 17 dB above each
% small scatterer (K ~ 6 times the rest combined) the exact law is even more
% concentrated than chi-square-4: the KS distance is 0.01 to the exact law
% and ~0.1 to chi-square-4.

%% [2026] What fluctuation costs in detection
Pfa = 1e-6; snr = 0:0.1:35;
pd0 = pd_swerling0(snr(1:5:end), Pfa); snr0 = snr(1:5:end);
pd1 = pd_swerling1(snr, Pfa); pd3 = pd_swerling3(snr, Pfa);
need = @(f) fzero(@(x) f(x) - 0.9, 15);
s0 = need(@(x) pd_swerling0(x, Pfa)); s1 = need(@(x) pd_swerling1(x, Pfa)); s3 = need(@(x) pd_swerling3(x, Pfa));
% Monte Carlo using the simulated RCS sequences of 5.5 and 5.6 as the target fluctuation
snr_mc = 10:4:30; n = 1e5; Y = -log(Pfa);
mc = zeros(2, numel(snr_mc)); seqs = {u12, u34};
for j = 1:2
    for i = 1:numel(snr_mc)
        uu = seqs{j}(randi(numel(seqs{j}), n, 1));
        z = abs(sqrt(db2lin(snr_mc(i))*uu).*exp(1i*2*pi*rand(n,1)) + sqrt(1/2)*(randn(n,1) + 1i*randn(n,1))).^2;
        mc(j, i) = mean(z > Y);
    end
end
fprintf('\n=== [2026] Fluctuation loss, single pulse, Pfa = 1e-6, Pd = 0.9 ===\n');
fprintf('Swerling 0 (constant): %.1f dB | Swerling I/II: %.1f dB (loss %.1f dB) | Swerling III/IV: %.1f dB (loss %.1f dB)\n', ...
        s0, s1, s1 - s0, s3, s3 - s0);
fig = figure('Name', '2026 Fluctuation loss');
plot(snr0, pd0, 'k', snr, pd1, 'r', snr, pd3, 'b', 'LineWidth', 1.5); hold on; grid on;
plot(snr_mc, mc(1,:), 'ro', snr_mc, mc(2,:), 'bs', 'MarkerSize', 7, 'LineWidth', 1.3);
yline(0.9, 'k:');
xlabel('Mean SNR [dB]'); ylabel('P_d'); title(sprintf('Single-pulse detection, P_{fa} = 10^{%d}', log10(Pfa)));
legend('Swerling 0 (constant RCS)', 'Swerling I/II', 'Swerling III/IV', ...
       'Monte Carlo, 5.5 RCS sequence', 'Monte Carlo, 5.6 RCS sequence', 'Location', 'southeast');
save_figure(fig, 's05_fluctuation_loss');
% A fluctuating target spends time in deep fades, so a single look needs much
% more mean SNR for high Pd: ~8 dB more for Swerling I at Pd = 0.9. A dominant
% scatterer (Swerling III) reduces that loss. Frequency agility or pulse-to-
% pulse decorrelation (Swerling II/IV) plus integration recover part of it.
