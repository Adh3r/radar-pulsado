%% DEMO_END_TO_END - Simulated pulsed-radar scene, from echoes to detections   [2026]
%  scene (radar equation + clutter + noise)  ->  data cube (fast x slow time)
%  -> Barker-13 matched filter  ->  fixed threshold vs CA-CFAR  ->  detections
%  -> Monte Carlo ROC versus theory
%
%  Builds on every lab: SNRs from the radar equation and the Lab 2 receiver
%  chain, the Barker-13 matched filter of Lab 4, thresholds and Pd from Lab 3,
%  and the cell-averaging CFAR of Lab 6. Slow-time clutter suppression (MTI,
%  Lab 7) is left as a documented extension point.

clear; close all;
run(fullfile(fileparts(mfilename('fullpath')), '..', 'startup.m'));

%% 1. Scene and data cube
sc  = scene_default();
out = simulate_scene(sc);
K   = numel(sc.code);
[nf, np] = size(out.rx);
r_km = out.range/1e3;
fprintf('\n=== Scene: %d fast-time bins x %d pulses, bin = %.0f m, Rmax = %.1f km ===\n', ...
        nf, np, out.range(2), out.range(end)/1e3);

%% 2. Slow-time processing
% Sections 3-6 run without MTI; section 6b repeats the chain with the single
% and double cancellers of Lab 7 in front of the matched filter.
cube = out.rx;

%% 3. Pulse compression: Barker-13 matched filter on every pulse
y = matched_filter(cube, sc.code, 'aligned') / sqrt(K);   % output noise power = 1
y(out.blind, :) = 0;

%% 4. Detection: fixed threshold (known noise power) vs CA-CFAR (Lab 6)
Pfa = 1e-6;
p1 = abs(y(:,1)).^2;            % single pulse
pN = sum(abs(y).^2, 2);         % non-coherent integration over the CPI
% fixed: assumes the noise power is known (= 1 here, by construction)
T1 = threshold_from_pfa(Pfa, 1, 1);
TN = threshold_from_pfa(Pfa, 1, np);
[d1, ~] = detect_peaks(p1, T1);
[dN, ~] = detect_peaks(pN, TN);
% CA-CFAR: noise estimated from M training cells, G guard cells per side
% (G = 3 keeps the compressed main lobe out of the training window)
cfar.M = 32; cfar.G = 3;
alpha1 = ca_cfar_alpha(Pfa, cfar.M, 1);
alphaN = ca_cfar_alpha(Pfa, cfar.M, np);      % cells integrate np pulses
[~, C1] = ca_cfar(p1, cfar.M, cfar.G, alpha1);
[~, CN] = ca_cfar(pN, cfar.M, cfar.G, alphaN);
[e1, ~] = detect_peaks(p1, C1);
[eN, ~] = detect_peaks(pN, CN);

%% 5. Score against the truth
binR = out.range(2);
inClutter = @(idx) out.range(idx) >= sc.clutter.r_min & out.range(idx) <= sc.clutter.r_max + K*binR;
hit = @(d, b) any(abs(d - b) <= 1);
status  = @(d, i) ternary(inClutter(out.bin(i)), 'masked', yesno(hit(d, out.bin(i))));
% inside the clutter belt a hit only counts as clean if no other detection lies within 300 m
others  = @(d, i) d(abs(d - out.bin(i)) > 1);
crowded = @(d, i) any(abs(others(d, i) - out.bin(i)) <= round(300/binR));
cstatus = @(d, i) ternary(inClutter(out.bin(i)) && hit(d, out.bin(i)) && crowded(d, i), 'yes*', yesno(hit(d, out.bin(i))));
fprintf('\n%-8s %6s %6s %9s %9s | %8s %8s | %8s %8s\n', 'target', 'R[km]', 'RCS', 'SNRchip', 'SNR_MF', ...
        'fix 1p', 'fix NCI', 'CFAR 1p', 'CFAR NCI');
for i = 1:size(sc.targets, 1)
    fprintf('%-8d %6.1f %6.2f %8.1f %9.1f | %8s %8s | %8s %8s\n', i, sc.targets(i,1)/1e3, sc.targets(i,2), ...
            out.snr_chip_dB(i), out.snr_chip_dB(i) + 10*log10(K), status(d1, i), status(dN, i), ...
            cstatus(e1, i), cstatus(eN, i));
end
fprintf('(* hit inside the clutter belt with other detections within 300 m: not separable from clutter)\n');
% false alarms: sidelobe (within one pulse length of a true echo), clutter, noise
classify = @(d) [nnz(arrayfun(@(x) min(abs(x - out.bin)) > 1 && min(abs(x - out.bin)) <= K, d) & ~inClutter(d)), ...
                 nnz(arrayfun(@(x) min(abs(x - out.bin)) > 1, d) & inClutter(d)), ...
                 nnz(arrayfun(@(x) min(abs(x - out.bin)) > K, d) & ~inClutter(d))];
fprintf('False alarms [sidelobe / clutter / noise] (noise-only expectation %.3f per profile):\n', Pfa*nf);
fprintf('  fixed  : 1 pulse %2d / %2d / %2d   NCI %2d / %2d / %2d\n', classify(d1), classify(dN));
fprintf('  CA-CFAR: 1 pulse %2d / %2d / %2d   NCI %2d / %2d / %2d   (M = %d, G = %d)\n', ...
        classify(e1), classify(eN), cfar.M, cfar.G);
% Reading: pulse compression lifts the weak 0.1 m^2 target at 12.3 km from
% ~0 dB to ~11 dB, where a single pulse detects it with Pd ~ 0.5 (a coin flip
% for any given realisation); integrating 32 pulses makes detection certain.
% Inside the clutter belt the fixed threshold fires on the clutter itself and
% the 0.01 m^2 drone cannot be separated from it: that is what MTI (Doppler
% separation from zero-Doppler clutter) and CFAR (locally estimated threshold)
% address. Integration has a side effect: the 32-pulse threshold sits only
% ~3 dB above the noise, while the -22 dB Barker range sidelobes of the
% 0.5 and 2 m^2 targets stay ~8-10 dB above it, so they are reported as false
% targets. Sidelobe control (mismatched filters, longer codes) and CFAR with
% guard cells are the standard answers.
%
% CA-CFAR (Lab 6), single pulse: estimating the noise locally removes almost
% all clutter false alarms (20 -> 1) and the sidelobe ones (7 -> 0), keeping
% every target outside the clutter. After integration it still cuts the
% sidelobe false targets (a strong peak inside the training window raises the
% threshold), but clutter false alarms GROW (3 -> 15): the clutter barely
% changes over the 3.2 ms CPI (sigma_f * CPI ~ 0.06), so an integrated clutter
% cell is not the Gamma(32) variable alpha was calibrated for, and the CFAR
% loses its constant-false-alarm property. Removing the zero-Doppler clutter
% before integration is the job of MTI (section 6b).

%% 6. Figures
clut = [sc.clutter.r_min, sc.clutter.r_max]/1e3;
fig = figure('Name', 'End-to-end demo', 'Position', [60 60 1100 850]);
panels = {20*log10(abs(out.rx(:,1)) + eps), 10*log10(p1 + eps), 10*log10(pN/np + eps)};
thr    = [NaN, 10*log10(T1), 10*log10(TN/np)];
dets   = {[], d1, dN};
ttl    = {'Received signal, pulse 1 (before compression): per-sample power [dB]', ...
          sprintf('Matched-filter output, pulse 1, fixed threshold P_{fa} = 10^{%d}', log10(Pfa)), ...
          sprintf('Non-coherent integration of %d pulses (mean power)', np)};
for k = 1:3
    ax = subplot(3,1,k); hold on; grid on; box on;
    patch([clut fliplr(clut)], [-20 -20 60 60], [1 0.85 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    plot(r_km, panels{k}, 'Color', [0.15 0.3 0.6]);
    if ~isnan(thr(k)), yline(thr(k), 'r--', 'threshold', 'LineWidth', 1.3); end
    if ~isempty(dets{k}), plot(r_km(dets{k}), panels{k}(dets{k}), 'ro', 'MarkerSize', 7, 'LineWidth', 1.4); end
    plot(sc.targets(:,1)/1e3, 55*ones(size(sc.targets,1),1), 'kv', 'MarkerFaceColor', [0.2 0.7 0.3], 'MarkerSize', 8);
    text(sc.targets(:,1)/1e3 + 0.15, 55*ones(size(sc.targets,1),1), ...
         compose('%.2f m^2', sc.targets(:,2)), 'FontSize', 8);
    ylim([-20 60]); xlim([0 r_km(end)]); ylabel('dB'); title(ttl{k});
    if k == 1, text(mean(clut), 48, 'clutter', 'HorizontalAlignment', 'center'); end
end
xlabel('Range [km]');
legend(ax, {'clutter belt', 'signal', 'threshold', 'detections', 'true targets'}, ...
       'Orientation', 'horizontal', 'Location', 'southoutside');
save_figure(fig, 'demo_overview');

fig = figure('Name', 'Fixed threshold vs CA-CFAR', 'Position', [60 60 1100 620]);
cpan = {p1, pN/np}; fixT = [T1, TN/np]; cT = {C1, CN/np}; cdet = {e1, eN};
cttl = {'Single pulse', sprintf('Non-coherent integration of %d pulses (mean power)', np)};
for k = 1:2
    ax = subplot(2,1,k); hold on; grid on; box on;
    patch([clut fliplr(clut)], [-20 -20 60 60], [1 0.85 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    plot(r_km, 10*log10(cpan{k} + eps), 'Color', [0.15 0.3 0.6]);
    yline(10*log10(fixT(k)), 'k--', 'LineWidth', 1.2);
    plot(r_km, 10*log10(cT{k}), 'r', 'LineWidth', 1.4);
    plot(r_km(cdet{k}), 10*log10(cpan{k}(cdet{k})), 'ro', 'MarkerSize', 7, 'LineWidth', 1.4);
    plot(sc.targets(:,1)/1e3, 55*ones(size(sc.targets,1),1), 'kv', 'MarkerFaceColor', [0.2 0.7 0.3], 'MarkerSize', 8);
    ylim([-20 60]); xlim([0 r_km(end)]); ylabel('dB');
    title(sprintf('%s: fixed threshold vs CA-CFAR (M = %d, G = %d), P_{fa} = 10^{%d}', ...
          cttl{k}, cfar.M, cfar.G, log10(Pfa)));
end
xlabel('Range [km]');
legend(ax, {'clutter belt', 'signal', 'fixed threshold', 'CA-CFAR threshold', 'CA-CFAR detections', 'true targets'}, ...
       'Orientation', 'horizontal', 'Location', 'southoutside');
save_figure(fig, 'demo_cfar');

fig = figure('Name', 'Range-pulse map');
imagesc(r_km, 1:np, 10*log10(abs(y.').^2 + eps)); axis xy; colorbar; caxis([-10 50]);
xlabel('Range [km]'); ylabel('Pulse #'); title('Matched-filter output over the CPI [dB]');
save_figure(fig, 'demo_range_pulse_map');

%% 6b. MTI before integration (Lab 7)
% Cancel the zero-Doppler clutter pulse to pulse, then matched filter,
% integrate the remaining pulses and detect with CA-CFAR. The canceller
% output is scaled by 1/sqrt(sum h^2) so that the noise power stays 1.
mti_cfg = {'no MTI', 0; 'single canceller', 1; 'double canceller', 2};
fprintf('\n=== MTI (Lab 7) + %d-pulse NCI + CA-CFAR (M = %d, G = %d) ===\n', np, cfar.M, cfar.G);
fprintf('%-17s | %5s %5s %5s %5s | %s\n', 'slow-time stage', 'T1', 'T2', 'T3', 'T4', ...
        'false alarms sidelobe / clutter / noise');
fig = figure('Name', 'MTI + CA-CFAR', 'Position', [60 60 1100 850]);
for c = 1:size(mti_cfg, 1)
    order = mti_cfg{c, 2};
    if order == 0
        z = out.rx;
    else
        [z, h] = mti_canceller(out.rx, order);
        z = z / sqrt(sum(h.^2));
    end
    yz = matched_filter(z, sc.code, 'aligned') / sqrt(K);
    yz(out.blind, :) = 0;
    nz = size(yz, 2);
    pz = sum(abs(yz).^2, 2);
    [~, Cz] = ca_cfar(pz, cfar.M, cfar.G, ca_cfar_alpha(Pfa, cfar.M, nz));
    dz = detect_peaks(pz, Cz);
    st = arrayfun(@(i) cstatus(dz, i), 1:size(sc.targets, 1), 'UniformOutput', false);
    fprintf('%-17s | %5s %5s %5s %5s | %2d / %2d / %2d\n', mti_cfg{c, 1}, st{:}, classify(dz));
    ax = subplot(3, 1, c); hold on; grid on; box on;
    patch([clut fliplr(clut)], [-20 -20 60 60], [1 0.85 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    plot(r_km, 10*log10(pz/nz + eps), 'Color', [0.15 0.3 0.6]);
    plot(r_km, 10*log10(Cz/nz), 'r', 'LineWidth', 1.4);
    plot(r_km(dz), 10*log10(pz(dz)/nz), 'ro', 'MarkerSize', 7, 'LineWidth', 1.4);
    plot(sc.targets(:,1)/1e3, 55*ones(size(sc.targets,1),1), 'kv', 'MarkerFaceColor', [0.2 0.7 0.3], 'MarkerSize', 8);
    ylim([-20 60]); xlim([0 r_km(end)]); ylabel('dB');
    title(sprintf('%s, then %d-pulse integration and CA-CFAR', mti_cfg{c, 1}, nz));
end
xlabel('Range [km]');
legend(ax, {'clutter belt', 'integrated power', 'CA-CFAR threshold', 'detections', 'true targets'}, ...
       'Orientation', 'horizontal', 'Location', 'southoutside');
save_figure(fig, 'demo_mti');
% Reading: the single canceller removes most of the clutter (CA ~ 38 dB for
% sigma_v = 1 m/s at 10 kHz PRF) and restores the CFAR behaviour after
% integration, and the drone (15 m/s, 0.03 PRF, -20 dB in the notch) survives.
% The double canceller removes the residual clutter completely, but its notch
% (-41 dB at the drone's Doppler) also removes the drone: the classic MTI
% trade-off between clutter rejection and slow-target visibility.

%% 7. Monte Carlo ROC: with and without pulse compression
snr_chip = out.snr_chip_dB(4);                   % weakest target (12.3 km, 0.1 m^2)
r = roc_monte_carlo(snr_chip, sc.code, 2e4, sc.seed + 1);
Pfa_axis = logspace(-5, 0, 60);
T = threshold_from_pfa(Pfa_axis, 1);
[pf_raw, pd_raw] = roc_empirical(r.h1_raw, r.h0_raw, T);
[pf_mf,  pd_mf ] = roc_empirical(r.h1_mf,  r.h0_mf,  T);
Pfa_th = logspace(-5, 0, 200);
fig = figure('Name', 'ROC');
semilogx(Pfa_th, pd_swerling0(snr_chip + r.gain_dB, Pfa_th), 'b-', 'LineWidth', 1.5); hold on;
semilogx(pf_mf, pd_mf, 'bo', 'MarkerSize', 5);
semilogx(Pfa_th, pd_swerling0(snr_chip, Pfa_th), 'r-', 'LineWidth', 1.5);
semilogx(pf_raw, pd_raw, 'rs', 'MarkerSize', 5);
grid on; xlabel('P_{fa}'); ylabel('P_d'); ylim([0 1]); xlim([1e-5 1]);
title(sprintf('ROC at %.1f dB per-sample SNR: Barker-13 compression adds %.1f dB', snr_chip, r.gain_dB));
legend('Compressed: theory', 'Compressed: Monte Carlo', 'Uncompressed: theory', ...
       'Uncompressed: Monte Carlo', 'Location', 'southeast');
save_figure(fig, 'demo_roc');
ok = pf_mf > 1e-4;
fprintf('ROC: max |Pd_MC - Pd_theory| (compressed, Pfa > 1e-4) = %.3f\n', ...
        max(abs(pd_mf(ok) - pd_swerling0(snr_chip + r.gain_dB, pf_mf(ok)))));

function s = yesno(b)
if b, s = 'yes'; else, s = 'no'; end
end

function v = ternary(c, a, b)
if c, v = a; else, v = b; end
end
