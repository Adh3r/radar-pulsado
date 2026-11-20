%% S03 - DETECTION THEORY: false alarms, Albersheim, exact Pd and Monte Carlo   (Lab 3)
%  Exercises 3.2-3.6 using src/detection. Extensions 2026: Albersheim is
%  checked against the exact Marcum-Q result and both against a Monte Carlo
%  simulation of an envelope detector; detection-range calculator combining
%  the Lab 2 radar equation with Albersheim.

clear; close all;
run(fullfile(fileparts(mfilename('fullpath')), '..', 'startup.m'));
rng(3);

%% 3.2 Mean time between false alarms vs threshold
ratio_dB = 0:0.1:15;
BW = 10.^(0:9);
Tfa = false_alarm_time(ratio_dB, BW);

fig = figure('Name', '3.2 Time between false alarms');
h = semilogy(ratio_dB, Tfa, 'LineWidth', 1.2); grid on;
col = turbo(numel(h)); for c = 1:numel(h), h(c).Color = col(c,:); end   % 10 distinct colours
xlabel('Threshold-to-noise ratio V_T^2/2\psi_0 [dB]'); ylabel('T_{fa} [s]');
title('Mean time between false alarms'); legend(compose('BW = 10^{%d} Hz', 0:9), 'Location', 'northwest');
for ref = [60 3600 86400 3.15e7], yline(ref, 'k:'); end
text(0.3*[1 1 1 1], [60 3600 86400 3.15e7]*1.8, {'1 min', '1 h', '1 day', '1 year'});
save_figure(fig, 's03_false_alarm_time');

%% 3.3 Effect of the threshold (BW = 1 MHz) and 3.4 effect of the bandwidth (12 dB)
[T3, P3] = false_alarm_time([8 10 12 14 15], 1e6);
fprintf('\n=== 3.3 Threshold (BW = 1 MHz) ===\n');
fprintf('%6s %10s %12s\n', 'dB', 'Pfa', 'Tfa [s]');
fprintf('%6.1f %10.2e %12.4g\n', [[8 10 12 14 15]; P3.'; T3.']);
T4 = false_alarm_time(12, [1 1e3 1e6 1e9]);
fprintf('\n=== 3.4 Bandwidth (threshold 12 dB) ===\n');
fprintf('BW = %8.0e Hz -> Tfa = %10.4g s\n', [[1 1e3 1e6 1e9]; T4]);
% Discussion
%  3.3 Tfa = exp(VT^2/2psi0)/BW is exponential in the LINEAR threshold:
%      +2 dB (12 -> 14 dB) multiplies Tfa by > 10^4. Raising VT also blocks
%      weak echoes, so VT is set from the tolerable Pfa and the SNR needed for
%      the desired Pd follows (3.5).
%  3.4 Tfa ~ 1/BW: each extra decade of bandwidth gives ten times more
%      independent noise samples per second, so ten times more false alarms.
%      Hence the IF filter is kept as narrow as the pulse allows (B ~ 1/tau).

%% 3.5-3.6 Detection curves with Albersheim
Pd  = 0.30:0.001:0.99;
Pfa = 10.^(-12:-3);
SNR_alb = albersheim_snr(Pd, Pfa);

fig = figure('Name', '3.6 Detection curves');
h = plot(SNR_alb, Pd, 'LineWidth', 1.2); grid on; ylim([0.3 1]);
col = turbo(numel(h)); for c = 1:numel(h), h(c).Color = col(c,:); end   % 10 distinct colours
xlabel('Single-pulse SNR [dB]'); ylabel('P_d');
title('Detection curves (Albersheim, one pulse, non-fluctuating target)');
legend(compose('P_{fa} = 10^{%d}', -12:-3), 'Location', 'southeast');
save_figure(fig, 's03_albersheim_curves');
fprintf('\n=== 3.6 Albersheim, Pfa = 1e-6 ===\n');
fprintf('Pd = %.2f -> SNR = %.2f dB\n', [[0.5 0.9 0.95 0.99]; albersheim_snr([0.5 0.9 0.95 0.99], 1e-6).']);
% Discussion: the curves are S-shaped and nearly parallel. Each decade of Pfa
% costs ~1 dB of SNR, but Pd gets expensive near 1 because B = ln(Pd/(1-Pd))
% diverges (Pd = 1 would need infinite SNR). Going from Pd 0.3 to 0.99 takes
% only ~4 dB, which is why a closing target goes from invisible to detected
% almost abruptly as its SNR rises with R^-4.

%% [2026] Exact Pd (Marcum Q) and Albersheim's error
Pd_grid  = [0.3 0.5 0.7 0.9 0.95 0.99];
Pfa_grid = 10.^(-12:-3);
SNR_grid = albersheim_snr(Pd_grid, Pfa_grid);
Pd_exact = pd_swerling0(SNR_grid, Pfa_grid);           % Pd that the Albersheim SNR actually gives
SNR_exact = zeros(size(SNR_grid));                     % exact SNR for each (Pd, Pfa)
for i = 1:numel(Pd_grid)
    for j = 1:numel(Pfa_grid)
        SNR_exact(i,j) = fzero(@(s) pd_swerling0(s, Pfa_grid(j)) - Pd_grid(i), SNR_grid(i,j));
    end
end
errA = SNR_grid - SNR_exact;
fprintf('\n=== [2026] Albersheim vs exact (Marcum Q) ===\n');
fprintf('max |error| = %.3f dB over Pd 0.3-0.99, Pfa 1e-12..1e-3\n', max(abs(errA(:))));
fprintf('Pd = 0.9, Pfa = 1e-6: Albersheim %.2f dB, exact %.2f dB\n', SNR_grid(4,7), SNR_exact(4,7));

fig = figure('Name', '2026 Albersheim error');
imagesc(log10(Pfa_grid), 1:numel(Pd_grid), errA); colorbar; axis xy;
set(gca, 'YTick', 1:numel(Pd_grid), 'YTickLabel', compose('%.2f', Pd_grid));
xlabel('log_{10} P_{fa}'); ylabel('P_d');
title('Albersheim SNR - exact SNR [dB]');
for i = 1:numel(Pd_grid)
    for j = 1:numel(Pfa_grid)
        text(log10(Pfa_grid(j)), i, sprintf('%.2f', errA(i,j)), 'HorizontalAlignment', 'center', 'FontSize', 8);
    end
end
save_figure(fig, 's03_albersheim_error');

%% [2026] Monte Carlo validation of Pfa and Pd (envelope / square-law detector)
n_noise = 2e6;   n_sig = 2e5;
noise = sqrt(1/2)*(randn(n_noise,1) + 1i*randn(n_noise,1));      % unit-power complex noise
Pfa_mc_target = 10.^(-5:-1);
T = threshold_from_pfa(Pfa_mc_target, 1);                          % square-law thresholds
Pfa_mc = mean(abs(noise).^2 > T, 1);

SNR_mc = 0:1:16;
Pfa_d  = 1e-4;  Td = threshold_from_pfa(Pfa_d, 1);
Pd_mc  = zeros(size(SNR_mc));
for i = 1:numel(SNR_mc)
    x = sqrt(db2lin(SNR_mc(i)))*exp(1i*2*pi*rand(n_sig,1)) + ...
        sqrt(1/2)*(randn(n_sig,1) + 1i*randn(n_sig,1));
    Pd_mc(i) = mean(abs(x).^2 > Td);
end
SNR_fine = 0:0.1:16;
fig = figure('Name', '2026 Monte Carlo validation', 'Position', [100 100 1000 400]);
subplot(1,2,1);
loglog(Pfa_mc_target, Pfa_mc, 'o', 'MarkerSize', 7, 'LineWidth', 1.3); hold on; grid on;
loglog([1e-5 1e-1], [1e-5 1e-1], 'k--');
xlabel('Design P_{fa} (threshold)'); ylabel('Measured P_{fa}');
title(sprintf('Rayleigh noise, %.0e samples', n_noise)); legend('Monte Carlo', 'Theory', 'Location', 'northwest');
subplot(1,2,2);
plot(SNR_fine, pd_swerling0(SNR_fine, Pfa_d), 'b-', 'LineWidth', 1.4); hold on; grid on;
plot(SNR_mc, Pd_mc, 'ro', 'LineWidth', 1.3);
Pd_alb_axis = 0.3:0.01:0.99;
plot(albersheim_snr(Pd_alb_axis, Pfa_d), Pd_alb_axis, 'k:', 'LineWidth', 1.4);
xlabel('SNR [dB]'); ylabel('P_d'); ylim([0 1]);
title(sprintf('P_d at P_{fa} = 10^{%d}', log10(Pfa_d)));
legend('Exact (Marcum Q)', sprintf('Monte Carlo (%.0e trials/pt)', n_sig), 'Albersheim', 'Location', 'southeast');
save_figure(fig, 's03_monte_carlo');

dev = abs(Pd_mc - pd_swerling0(SNR_mc, Pfa_d));
fprintf('\n=== [2026] Monte Carlo ===\n');
fprintf('max |Pfa_mc/Pfa - 1| = %.3f   (1e-5 point: %.0f expected crossings)\n', ...
        max(abs(Pfa_mc./Pfa_mc_target - 1)), 1e-5*n_noise);
fprintf('max |Pd_mc - Pd_exact| = %.4f (binomial std <= %.4f)\n', max(dev), sqrt(0.25/n_sig));

%% [2026] Detection range: radar equation (Lab 2) + Albersheim (Lab 3)
% Lab 2 scenario: 1.5 MW, 1.1 GHz, 50 MHz, 1 m^2, 30 dBi (assumed), both receiver chains.
Pt = 1.5e6; G_ant = 30; fc = 1.1e9; B = 50e6; sigma = 1;
NFA = cascade_noise_figure([37 -2.7 -1 -6 -1 11], [1.2 2.7 1 5.5 1 1.1]);   % LNA first
NFB = cascade_noise_figure([-2.7 37 -1 -6 -1 11], [2.7 1.2 1 5.5 1 1.1]);   % cable first
Pd  = [0.5 0.8 0.9 0.95 0.99];
Pfa = 10.^(-8:-4);
SNRreq = albersheim_snr(Pd, Pfa);                                 % Pd x Pfa
Rdet = max_detection_range(Pt, G_ant, fc, B, NFA, 0, sigma, SNRreq);

fig = figure('Name', '2026 Detection range');
plot(SNRreq(:, 3), Rdet(:, 3)/1e3, 'o-', 'LineWidth', 1.5); hold on; grid on;
for j = [1 5]
    plot(SNRreq(:, j), Rdet(:, j)/1e3, '.--');
end
text(SNRreq(:,3) + 0.1, Rdet(:,3)/1e3 + 1.5, compose('P_d = %.2f', Pd));
xlabel('Required single-pulse SNR [dB]'); ylabel('Detection range [km]');
title(sprintf('1 m^2 target, chain A, G = %d dBi: detection range vs (P_d, P_{fa})', G_ant));
legend('P_{fa} = 10^{-6}', 'P_{fa} = 10^{-8}', 'P_{fa} = 10^{-4}', 'Location', 'northeast');
save_figure(fig, 's03_detection_range');
fprintf('\n=== [2026] Detection range for a 1 m^2 target ===\n');
fprintf('Pd = 0.9, Pfa = 1e-6: SNR %.1f dB -> R = %.1f km (chain A) / %.1f km (chain B)\n', ...
        SNRreq(3,3), Rdet(3,3)/1e3, ...
        max_detection_range(Pt, G_ant, fc, B, NFB, 0, sigma, SNRreq(3,3))/1e3);
