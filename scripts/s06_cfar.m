%% S06 - CFAR: cell-averaging constant false alarm rate detector   (Lab 6)
%  Exercises 6.1-6.4 using src/detection (ca_cfar, ca_cfar_alpha, ca_cfar_pfa).
%  Extensions 2026: large Monte Carlo with confidence intervals, the
%  linear-detector pitfall, CFAR loss for Swerling I targets, and the two
%  classic CA-CFAR failure modes (clutter edge, multiple targets).
%
%  A fixed threshold only works if the noise power is known. CA-CFAR
%  estimates it locally, from M training cells around the cell under test
%  (CUT), and scales that estimate by alpha so that Pfa stays constant
%  whatever the noise level.

clear; close all;
run(fullfile(fileparts(mfilename('fullpath')), '..', 'startup.m'));
rng(6);
k = radar_constants();

%% 6.1 Noise statistics before and after the detector
L  = 10000;
T0 = 300; B = 1e6;                         % lab values
sigma2 = k.kB*T0*B;                        % noise power of EACH component [W]
xI = sqrt(sigma2)*randn(L,1);
xQ = sqrt(sigma2)*randn(L,1);
y_sq  = xI.^2 + xQ.^2;                     % square-law detector (power)
y_lin = sqrt(y_sq);                        % linear detector (envelope)

fprintf('\n=== 6.1 Detector statistics (kT0B = %.3e W) ===\n', sigma2);
fprintf('power of xI: %.3e W, xQ: %.3e W\n', mean(xI.^2), mean(xQ.^2));
fprintf('square law: mean %.3e (theory 2kT0B = %.3e), std/mean %.3f (exponential: 1)\n', ...
        mean(y_sq), 2*sigma2, std(y_sq)/mean(y_sq));
fprintf('linear:     mean %.3e (theory sqrt(pi/2) sqrt(kT0B) = %.3e)\n', mean(y_lin), sqrt(pi/2*sigma2));

s = sqrt(sigma2);
fig = figure('Name', '6.1 Detector statistics', 'Position', [100 100 1150 360]);
subplot(1,3,1);
histogram(xI/s, 60, 'Normalization', 'pdf', 'FaceColor', [0.6 0.7 0.9]); hold on; grid on;
u = linspace(-4.5, 4.5, 300); plot(u, exp(-u.^2/2)/sqrt(2*pi), 'r', 'LineWidth', 1.5);
xlabel('x_I / \sigma'); ylabel('pdf'); title('Input (I channel): Gaussian');
subplot(1,3,2);
histogram(y_sq/(2*sigma2), 60, 'Normalization', 'pdf', 'FaceColor', [0.6 0.7 0.9]); hold on; grid on;
u = linspace(0, 8, 300); plot(u, exp(-u), 'r', 'LineWidth', 1.5); xlim([0 8]);
xlabel('y / 2\sigma^2'); title('Square-law x_I^2 + x_Q^2: exponential');
subplot(1,3,3);
histogram(y_lin/s, 60, 'Normalization', 'pdf', 'FaceColor', [0.6 0.7 0.9]); hold on; grid on;
u = linspace(0, 4.5, 300); plot(u, u.*exp(-u.^2/2), 'r', 'LineWidth', 1.5);
xlabel('r / \sigma'); title('Linear (x_I^2 + x_Q^2)^{1/2}: Rayleigh');
legend('Simulated', 'Theory');
save_figure(fig, 's06_detector_statistics');
% Discussion: I and Q are independent zero-mean Gaussians. The square-law
% output is a chi-square with 2 degrees of freedom, i.e. an exponential
% (std = mean), and the linear-detector output is Rayleigh. Eq. 6.1 is
% derived for the exponential (square-law) samples. (The statement calls
% x_I^2 + x_Q^2 the envelope; strictly it is the power, the envelope being
% its square root.)

%% 6.2 Scaling factor for M = 40 and Pfa = 0.01
M = 40; Pfa0 = 0.01;
alpha0 = ca_cfar_alpha(Pfa0, M);
fprintf('\n=== 6.2 Scaling factor ===\n');
fprintf('M = %d, Pfa = %g -> alpha = M (Pfa^(-1/M) - 1) = %.4f\n', M, Pfa0, alpha0);
fprintf('(known noise power would need alpha = -ln(Pfa) = %.4f; the excess is the price of estimating it)\n', -log(Pfa0));

%% 6.3 CA-CFAR along a range profile
G = 2;                                         % guard cells per side (the lab leaves G free)
[det, T] = ca_cfar(y_sq, M, G, alpha0);
valid = ~isnan(T);
n_fa = nnz(det);
fprintf('\n=== 6.3 CA-CFAR (M = %d, G = %d per side), %d valid cells ===\n', M, G, nnz(valid));
fprintf('false alarms: %d -> Pfa = %.4f (design %.2f)\n', n_fa, n_fa/nnz(valid), Pfa0);

idx = 1:400;
fig = figure('Name', '6.3 CA-CFAR', 'Position', [100 100 1100 380]);
semilogy(idx, y_sq(idx)/(2*sigma2), 'Color', [0.2 0.35 0.65]); hold on; grid on;
semilogy(idx, T(idx)/(2*sigma2), 'r', 'LineWidth', 1.5);
fa = find(det(idx));
semilogy(fa, y_sq(fa)/(2*sigma2), 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 7);
xlabel('Range cell'); ylabel('Power / noise power');
title(sprintf('CUT and CA-CFAR threshold (M = %d, \\alpha = %.2f): %d false alarms in %d cells, P_{fa} = %.4f', ...
      M, alpha0, n_fa, nnz(valid), n_fa/nnz(valid)));
legend('CUT (square-law samples)', 'Adaptive threshold', 'False alarms', 'Location', 'southeast');
save_figure(fig, 's06_cfar_threshold');
% The threshold follows the local noise estimate, so it fluctuates; only the
% scaling factor fixes Pfa.

%% 6.4 Check the routine for Pfa = 0.1, 0.01, 0.001 and 0.0001
% Window-matrix simulation suggested in the lab: each column is one
% independent window (training + guard + CUT), the threshold is the product
% of a coefficient row vector with that matrix.
Pfa_set = [0.1 0.01 1e-3 1e-4];
W  = M + 2*G + 1;                               % window length
w  = [ones(1,M/2), zeros(1,G), 0, zeros(1,G), ones(1,M/2)] / M;   % coefficient vector
cut = M/2 + G + 1;
mc = @(n) -log(rand(W, n));                     % exponential samples, unit mean
fprintf('\n=== 6.4 Window-matrix Monte Carlo ===\n');
fprintf('%9s %10s %10s %10s %12s %22s\n', 'Pfa', 'trials', 'expected', 'counted', 'Pfa_est', '95% interval');
res = zeros(2, numel(Pfa_set));
for pass = 1:2
    n = [1e4, 2e6]; n = n(pass);                % lab size, then [2026] large run
    X = mc(n);
    noise = w * X;                              % 1 x n training means
    for i = 1:numel(Pfa_set)
        c = nnz(X(cut, :) > ca_cfar_alpha(Pfa_set(i), M) * noise);
        res(pass, i) = c/n;
        [lo, hi] = binomial_interval(c, n);
        fprintf('%9.0e %10.0e %10.1f %10d %12.2e   [%8.2e, %8.2e]\n', Pfa_set(i), n, Pfa_set(i)*n, c, c/n, lo, hi);
    end
    if pass == 1, fprintf('--- [2026] same test with 2e6 trials ---\n'); end
end
% Discussion (6.4): with the lab size of 10^4 trials the estimates are only
% indicative: 0.1 and 0.01 land within ~20 % (and can fall just outside the
% 95 % interval by chance), 10^-3 rests on ~10 events and 10^-4 on ~1 event,
% so the latter is dominated by chance (0, 1 or 2 counts). The 2e6-trial run
% puts every design value inside its 95 % interval. Before blaming the CFAR,
% suggested next steps: (1) increase the trials to >= 100/Pfa (relative error
% ~ 1/sqrt(n Pfa)) or use importance sampling, (2) check that the samples are
% square-law (eq. 6.1 is not valid for a linear detector, see below), and
% (3) check that training cells are independent (oversampled or filtered data
% reduces the effective M).

%% [2026] Pitfall: the same alpha with a linear detector
Xl = sqrt(mc(2e6));                                       % Rayleigh envelopes
fprintf('\n=== [2026] Linear detector with the square-law alpha ===\n');
for i = 1:numel(Pfa_set)
    fprintf('design %.0e -> measured %.2e\n', Pfa_set(i), ...
            mean(Xl(cut,:) > ca_cfar_alpha(Pfa_set(i), M) * (w*Xl)));
end
% Applying the square-law alpha to envelope samples gives a completely
% different (much lower) Pfa: the threshold must match the detector law.

fig = figure('Name', '6.4 Pfa check');
loglog(Pfa_set, res(1,:), 'rs', 'MarkerSize', 9, 'LineWidth', 1.5); hold on; grid on;
loglog(Pfa_set, res(2,:), 'bo', 'MarkerSize', 7, 'LineWidth', 1.5);
loglog([5e-5 0.2], [5e-5 0.2], 'k--');
xlabel('Design P_{fa}'); ylabel('Measured P_{fa}');
title(sprintf('CA-CFAR, M = %d: measured vs design P_{fa}', M));
legend('10^4 trials (lab)', '2\cdot10^6 trials [2026]', 'Ideal', 'Location', 'northwest');
save_figure(fig, 's06_pfa_check');

%% [2026] CFAR loss for a Swerling I target
Pfa_d = 1e-6; snr = 0:0.25:35; Ms = [8 16 40 Inf];
fig = figure('Name', '2026 CFAR loss');
fprintf('\n=== [2026] CFAR loss, Swerling I, Pfa = 1e-6, Pd = 0.5 ===\n');
snr50_fixed = fzero(@(x) pd_swerling1(x, Pfa_d) - 0.5, 10);
for Mi = Ms
    plot(snr, pd_swerling1(snr, Pfa_d, Mi), 'LineWidth', 1.5); hold on;
    if ~isinf(Mi)
        s50 = fzero(@(x) pd_swerling1(x, Pfa_d, Mi) - 0.5, 15);
        fprintf('M = %2d: SNR(Pd=0.5) = %.2f dB, CFAR loss = %.2f dB\n', Mi, s50, s50 - snr50_fixed);
    end
end
% Monte Carlo through ca_cfar's window for M = 40
snr_mc = 10:5:30; n = 2e5; pd_mc = zeros(size(snr_mc)); a40 = ca_cfar_alpha(Pfa_d, 40);
for i = 1:numel(snr_mc)
    X = mc(n); X(cut,:) = X(cut,:) * (1 + db2lin(snr_mc(i)));   % Swerling I: exponential with mean 1+SNR
    pd_mc(i) = mean(X(cut,:) > a40 * (w*X));
end
plot(snr_mc, pd_mc, 'ko', 'MarkerSize', 7, 'LineWidth', 1.4);
grid on; xlabel('SNR [dB]'); ylabel('P_d');
title(sprintf('Swerling I, P_{fa} = 10^{%d}: CA-CFAR vs known noise power', log10(Pfa_d)));
legend('M = 8', 'M = 16', 'M = 40', 'Fixed threshold (M \rightarrow \infty)', 'Monte Carlo, M = 40', ...
       'Location', 'southeast');
save_figure(fig, 's06_cfar_loss');
fprintf('Monte Carlo M = 40: max |Pd_MC - Pd_theory| = %.4f\n', max(abs(pd_mc - pd_swerling1(snr_mc, Pfa_d, 40))));

%% [2026] Where CA-CFAR fails: target masking and clutter edges
% (a) one realisation: a strong target masks a weak neighbour, and a target
%     just before a clutter edge is masked by the clutter in its window
n_cells = 400; edge = 260; Pfa_e = 1e-3;
a = ca_cfar_alpha(Pfa_e, M);
p = -log(rand(n_cells, 1));                                % unit-mean noise
p(edge:end) = p(edge:end) * db2lin(20);                    % +20 dB clutter from the edge on
tg = [100 106 250];  snr_t = [30 13 13];                   % strong, weak neighbour, weak near edge
p(tg) = (1 + db2lin(snr_t(:))) .* (-log(rand(3, 1)));      % Swerling I samples
[d2, T2] = ca_cfar(p, M, G, a);
% (b) Monte Carlo: local Pfa of every cell around the clutter edge
n_mc = 2e4;
P = -log(rand(n_cells, n_mc));  P(edge:end, :) = P(edge:end, :) * db2lin(20);
pfa_local = mean(ca_cfar(P, M, G, a), 2);
cells = (M/2 + G + 1):(n_cells - M/2 - G);
fprintf('\n=== [2026] Failure modes (M = %d, Pfa = %.0e) ===\n', M, Pfa_e);
fprintf('targets %s (SNR %s dB): detected %s\n', mat2str(tg), mat2str(snr_t), mat2str(d2(tg)'));
fprintf('noise est. at the weak neighbour: %.1f x noise (strong target in its window)\n', T2(106)/a);
fprintf('local Pfa: %.1e far from the edge, max %.1e just inside the clutter, min %.1e just before it\n', ...
        mean(pfa_local(cells(1):edge-M/2-G-1)), max(pfa_local(cells)), min(pfa_local(edge-M/2-G:edge-1)));

fig = figure('Name', '2026 CFAR failure modes', 'Position', [100 100 1100 700]);
subplot(2,1,1);
semilogy(p, 'Color', [0.2 0.35 0.65]); hold on; grid on;
semilogy(T2, 'r', 'LineWidth', 1.5);
semilogy(find(d2), p(d2), 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 7);
semilogy(tg, p(tg), 'gv', 'MarkerFaceColor', 'g', 'MarkerSize', 9);
text(tg + 4, p(tg)*2, compose('%d dB', snr_t));
xline(edge, 'k--', 'clutter edge (+20 dB)', 'LabelVerticalAlignment', 'bottom');
ylabel('Power / noise power'); title('(a) Masking: strong neighbour and clutter edge raise the threshold');
legend('CUT', 'CA-CFAR threshold', 'Detections', 'True targets', 'Location', 'northwest');
subplot(2,1,2);
semilogy(cells, max(pfa_local(cells), 1/n_mc), 'b', 'LineWidth', 1.4); hold on; grid on;
yline(Pfa_e, 'k--', 'design P_{fa}'); xline(edge, 'k--');
xlabel('Range cell'); ylabel('Local P_{fa}'); ylim([1/n_mc 1]);
title(sprintf('(b) Local false-alarm rate around the edge (%.0e trials per cell)', n_mc));
save_figure(fig, 's06_cfar_failure_modes');
% Discussion: a strong target inside the training window inflates the noise
% estimate and masks a weak neighbour; a target just before a clutter edge
% sees clutter in half of its window (threshold too high: masked), and cells
% just inside the clutter see noise in half of theirs (threshold ~alpha/2
% times the clutter level: Pfa jumps by more than an order of magnitude).
% GO-CFAR (greatest-of the two half windows) handles edges; SO-CFAR and
% OS-CFAR (ordered statistic) handle multiple targets.

%% Local function
function [lo, hi] = binomial_interval(c, n)
% 95 % Wilson score interval for c successes in n trials.
z = 1.96; ph = c/n;
centre = (ph + z^2/(2*n)) / (1 + z^2/n);
half = z*sqrt(ph*(1 - ph)/n + z^2/(4*n^2)) / (1 + z^2/n);
lo = max(0, centre - half); hi = centre + half;
end
