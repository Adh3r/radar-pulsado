%% S07 - MTI: single and double delay-line cancellers   (Lab 7)
%  Exercises 7.1-7.4 using src/processing (mti_canceller, mti_response).
%  Extensions 2026: clutter attenuation versus clutter spread (eq. 7.1),
%  simulated and closed form, and blind speeds of the demo radar.
%
%  Stationary clutter returns the same phase on every pulse, a moving target
%  does not. Subtracting consecutive pulses (a high-pass filter in slow time)
%  removes what does not change: the clutter.

clear; close all;
run(fullfile(fileparts(mfilename('fullpath')), '..', 'startup.m'));
rng(7);
L = 10001;                                     % lab length
fnorm = @(N) (0:N-1)/N;                        % Doppler axis normalised to the PRF (f PRI)

%% 7.1 White Gaussian noise: flat spectrum, Gaussian histogram
x = randn(L, 1);                               % zero mean, 1 W
X2 = abs(fft(x)).^2 / L;                       % periodogram (mean = noise power)
fprintf('\n=== 7.1 White noise ===\n');
fprintf('power %.4f W, mean %.4f; periodogram mean %.4f, std %.4f (exponential-like: std ~ mean)\n', ...
        mean(x.^2), mean(x), mean(X2), std(X2));
fig = figure('Name', '7.1 White noise', 'Position', [100 100 1000 360]);
subplot(1,2,1);
plot(fnorm(L), 10*log10(X2), 'Color', [0.6 0.7 0.9]); hold on; grid on;
plot(fnorm(L), 10*log10(movmean(X2, 201)), 'b', 'LineWidth', 1.5); yline(0, 'r--', 'LineWidth', 1.3);
xlabel('f \cdot PRI'); ylabel('dB'); ylim([-40 15]);
title('Spectrum: flat at the noise power (0 dB)'); legend('periodogram', 'moving average', '1 W', 'Location', 'southwest');
subplot(1,2,2);
histogram(x, 60, 'Normalization', 'pdf', 'FaceColor', [0.6 0.7 0.9]); hold on; grid on;
u = linspace(-4.5, 4.5, 300); plot(u, exp(-u.^2/2)/sqrt(2*pi), 'r', 'LineWidth', 1.5);
xlabel('x'); ylabel('pdf'); title('Histogram: Gaussian N(0, 1)'); legend('simulated', 'theory');
save_figure(fig, 's07_white_noise');
% A single periodogram fluctuates around the true level (each bin is
% exponential with std = mean): averaging is needed to see the flat PSD.

%% 7.2 Single canceller on the noise vector
% (the statement says "double canceller" here; the title and 7.3 refer to
%  the single one, which is what is simulated)
[y1, h1] = mti_canceller(x, 1);
Y2 = abs(fft(y1)).^2 / numel(y1);
fprintf('\n=== 7.2 Single canceller ===\n');
fprintf('output power %.3f W (theory sum(h^2) = %d W)\n', mean(y1.^2), sum(h1.^2));
fig = figure('Name', '7.2 Single canceller', 'Position', [100 100 1000 360]);
subplot(1,2,1);
n = 1:60;
plot(n, x(n+1), 'b.-', n, y1(n), 'r.-'); grid on;
xlabel('Pulse #'); ylabel('Amplitude'); legend('input x[n]', 'output x[n] - x[n-1]');
title('Time domain: output power doubles');
subplot(1,2,2);
f2 = fnorm(numel(y1));
plot(fnorm(L), 10*log10(movmean(X2, 201)), 'b', 'LineWidth', 1.3); hold on; grid on;
plot(f2, 10*log10(movmean(Y2, 201)), 'r', 'LineWidth', 1.3);
plot(f2, 20*log10(mti_response(f2, 1, 1) + eps), 'k:', 'LineWidth', 1.4);
xlabel('f \cdot PRI'); ylabel('dB'); ylim([-30 10]);
legend('input (smoothed)', 'output (smoothed)', '|H_{SC}|^2 theory', 'Location', 'south');
title('Spectrum: notch at 0 and at the PRF');
save_figure(fig, 's07_single_canceller');
% Yes: in time the output looks like noise of twice the power (sum h^2 = 2);
% in frequency it is the input shaped by 4 sin^2(pi f PRI), with nulls at
% f = 0 and f = PRF (where the stationary clutter would be).

%% 7.3-7.4 Average frequency response over 1000 runs vs eqs. 7.2-7.3
n_runs = 1000;
fprintf('\n=== 7.3-7.4 Averaged responses (%d runs) ===\n', n_runs);
fig = figure('Name', '7.3-7.4 Canceller responses', 'Position', [100 100 1000 380]);
cols = {[0.85 0.33 0.1], [0 0.45 0.74]};
for order = 1:2
    N = L - order;
    acc_in = zeros(N, 1); acc_out = zeros(N, 1);
    for r = 1:n_runs
        xr = randn(L, 1);
        yr = mti_canceller(xr, order);
        acc_in  = acc_in  + abs(fft(xr(order+1:end))).^2;      % input aligned with the output
        acc_out = acc_out + abs(fft(yr)).^2;
    end
    Hest = acc_out ./ acc_in;                                  % averaged |H|^2
    Hest = Hest / max(Hest);                                   % normalised to 1 (as asked)
    fN = fnorm(N);
    Hth = (mti_response(fN, 1, order) / 2^order).^2;           % normalised |H|^2
    err = max(abs(Hest(:) - Hth(:)));
    fprintf('order %d: max |simulated - theory| (normalised |H|^2) = %.4f\n', order, err);
    subplot(1,2,1); plot(fN, Hest, 'Color', cols{order}, 'LineWidth', 1.1); hold on;
    plot(fN, Hth, 'k:', 'LineWidth', 1.4);
    subplot(1,2,2); plot(fN, 10*log10(Hest), 'Color', cols{order}, 'LineWidth', 1.1); hold on;
    plot(fN, 10*log10(Hth + eps), 'k:', 'LineWidth', 1.4);
end
subplot(1,2,1); grid on; xlabel('f \cdot PRI'); ylabel('Normalised |H|^2');
title('Linear'); legend('single (sim)', 'sin^2(\pi f PRI)', 'double (sim)', 'sin^4(\pi f PRI)', 'Location', 'south');
subplot(1,2,2); grid on; xlabel('f \cdot PRI'); ylabel('dB'); ylim([-60 3]);
title('dB: the double canceller has a wider, deeper notch');
save_figure(fig, 's07_canceller_response');
% Discussion: the averaged estimates follow sin^2 and sin^4 (eqs. 7.2 and 7.3
% squared, normalised). The double canceller widens the stop band around
% zero Doppler, which rejects spread clutter better, at the price of
% attenuating slow targets more and of a narrower pass band.

%% [2026] Clutter attenuation vs clutter spread (Gaussian spectrum, eq. 7.1)
PRI = 1; sf = logspace(-3, log10(0.2), 25);          % sigma_f * PRI
n_cells = 4000; n_pulses = 3 + 2;
CA_sim = zeros(2, numel(sf));
for i = 1:numel(sf)
    c = clutter_slow_time(n_cells, n_pulses, sf(i), PRI, 1);
    for order = 1:2
        CA_sim(order, i) = 10*log10(mean(abs(c(:)).^2) / mean(abs(reshape(mti_canceller(c, order), [], 1)).^2));
    end
end
CA_th = [mti_clutter_attenuation(sf, PRI, 1); mti_clutter_attenuation(sf, PRI, 2)];
[~, I1] = mti_clutter_attenuation(0.01, 1, 1); [~, I2] = mti_clutter_attenuation(0.01, 1, 2);
fprintf('\n=== [2026] Clutter attenuation, Gaussian clutter ===\n');
fprintf('sigma_f PRI = 0.01: CA single %.1f dB, double %.1f dB; improvement factor %.1f / %.1f dB\n', ...
        mti_clutter_attenuation(0.01, 1, 1), mti_clutter_attenuation(0.01, 1, 2), I1, I2);
fprintf('max |CA_sim - CA_theory| = %.2f dB\n', max(abs(CA_sim(:) - CA_th(:))));

fig = figure('Name', '2026 Clutter attenuation', 'Position', [100 100 1100 380]);
subplot(1,2,1);
f = linspace(-1, 1, 2001); sfe = 0.05;
Sc = zeros(size(f)); for kk = -2:2, Sc = Sc + exp(-(f - kk).^2/(2*sfe^2)); end
plot(f, Sc/max(Sc), 'Color', [0.6 0.35 0.1], 'LineWidth', 1.5); hold on; grid on;
plot(f, (mti_response(f, 1, 1)/2).^2, 'Color', cols{1}, 'LineWidth', 1.3);
plot(f, (mti_response(f, 1, 2)/4).^2, 'Color', cols{2}, 'LineWidth', 1.3);
xlabel('f \cdot PRI'); ylabel('Normalised'); title(sprintf('Clutter spectrum (\\sigma_f PRI = %.2f) and MTI responses', sfe));
legend('clutter (periodic in PRF)', 'single canceller', 'double canceller', 'Location', 'south');
subplot(1,2,2);
semilogx(sf, CA_th(1,:), '-', 'Color', cols{1}, 'LineWidth', 1.5); hold on; grid on;
semilogx(sf, CA_th(2,:), '-', 'Color', cols{2}, 'LineWidth', 1.5);
semilogx(sf, CA_sim(1,:), 'o', 'Color', cols{1}); semilogx(sf, CA_sim(2,:), 's', 'Color', cols{2});
xlabel('\sigma_f \cdot PRI'); ylabel('Clutter attenuation [dB]');
title('CA = P_{c,in} / P_{c,out}: closed form (lines), simulation (markers)');
legend('single', 'double', 'Location', 'northeast');
save_figure(fig, 's07_clutter_attenuation');
% Every decade of narrower clutter spectrum buys 20 dB (single) or 40 dB
% (double) of attenuation: the cancellers are 1st and 2nd order high-pass
% filters. Real clutter is never perfectly still, so CA is finite.

%% [2026] Blind speeds of the demo radar
sc = scene_default(); lambda = 3e8/sc.fc; PRF = 1/sc.PRI;
v = linspace(-1.2, 1.2, 4001) * lambda*PRF/2 * 2;           % +-2.4 blind speeds
fdv = doppler_frequency(v, sc.fc);
fig = figure('Name', '2026 Blind speeds');
plot(v, 20*log10(mti_response(fdv, sc.PRI, 1)/2 + eps), 'Color', cols{1}, 'LineWidth', 1.3); hold on; grid on;
plot(v, 20*log10(mti_response(fdv, sc.PRI, 2)/4 + eps), 'Color', cols{2}, 'LineWidth', 1.3);
for t = 1:size(sc.targets, 1)
    xline(sc.targets(t, 3), 'k--', sprintf('T%d', t), 'LabelVerticalAlignment', 'bottom');
end
ylim([-60 3]); xlabel('Radial velocity [m/s]'); ylabel('Normalised power gain [dB]');
title(sprintf('Demo radar (f_c = %.0f GHz, PRF = %.0f kHz): blind speeds every %.0f m/s', ...
      sc.fc/1e9, PRF/1e3, lambda*PRF/2));
legend('single canceller', 'double canceller', 'Location', 'southwest');
save_figure(fig, 's07_blind_speeds');
fprintf('\n=== [2026] Blind speeds: v_b = lambda PRF / 2 = %.0f m/s ===\n', lambda*PRF/2);
for t = 1:size(sc.targets, 1)
    fd_t = doppler_frequency(sc.targets(t,3), sc.fc);
    fprintf('target %d (%4.0f m/s, fd %5.0f Hz): gain single %6.1f dB, double %6.1f dB (normalised)\n', t, ...
            sc.targets(t,3), fd_t, 20*log10(mti_response(fd_t, sc.PRI, 1)/2), 20*log10(mti_response(fd_t, sc.PRI, 2)/4));
end
% The drone (15 m/s, 300 Hz = 0.03 PRF) sits deep in the notch: -20 dB with
% the single and -41 dB with the double canceller. Slow targets are the
% price of MTI; staggered PRFs move the blind speeds.
