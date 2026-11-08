%% S01 - PULSED RADAR: pulse burst, range ambiguity and spectrum   (Lab 1)
%  Exercises 1.2-1.4 using src/waveform. Extension 2026: the FFT spectra are
%  checked against the closed-form sinc x Dirichlet envelope.
%
%  A pulsed radar measures the round-trip time of short RF bursts: R = c t / 2.
%  tau, PRI and fc fix what it can and cannot measure.

clear; close all;
run(fullfile(fileparts(mfilename('fullpath')), '..', 'startup.m'));
k = radar_constants();

%% 1.2 Pulse burst with the lab parameters
tau = 1e-6; PRI = 5e-6; fc = 5e6; Pt_dBW = 60; L_dB = 70; Np = 5; fs = 100e6;
[m, v, t, env] = pulsed_radar(tau, PRI, fc, Pt_dBW, L_dB, Np, fs, 0);

fprintf('\n=== 1.2 Radar figures ===\n');
fprintf('PRF %.0f kHz | blind range %.0f m | unambiguous range %.0f m | resolution %.0f m\n', ...
        m.PRF/1e3, m.Rmin, m.Rmax, m.dR);
fprintf('duty %.0f %% | mean Rx power %.2f dBm | pulse energy %.2e J | peak voltage %.3f V\n', ...
        m.duty, m.Pavg_dBm, m.Epulse, m.Vpeak);
% Mean power = peak power x duty: -10 dBW = +20 dBm peak, minus 10log10(5) = 7 dB.

fig = figure('Name', '1.2 Pulse burst');
plot(t*1e6, v, 'Color', [0 0.35 0.75]); hold on; grid on;
plot(t*1e6, m.Vpeak*env, 'r--', t*1e6, -m.Vpeak*env, 'r--', 'LineWidth', 1.1);
xlabel('Time [\mus]'); ylabel('Voltage [V]');
title(sprintf('%d pulses: \\tau = %.0f \\mus, PRI = %.0f \\mus, f_c = %.0f MHz (5 carrier cycles per pulse)', ...
      Np, tau*1e6, PRI*1e6, fc/1e6));
legend('RF signal', 'Envelope', 'Location', 'best');
save_figure(fig, 's01_pulse_burst');

%% 1.3 Two targets and range ambiguity
R  = [450 1852];             % target 1 and target 2 (1 NM)
Lt = [3 6];                  % attenuation relative to the transmitted pulse [dB]
td = 2*R/k.c;
Rapp = k.c*mod(td, PRI)/2;   % range the radar reports
fprintf('\n=== 1.3 Two targets (unambiguous range %.0f m) ===\n', m.Rmax);
for i = 1:2
    fprintf('Target %d: R = %4.0f m, delay %.2f us = %.2f PRI -> measured %4.0f m\n', ...
            i, R(i), td(i)*1e6, td(i)/PRI, Rapp(i));
end

% Timeline of transmitted pulses and labelled echoes over 8 PRIs
Nvis = 8; tv = (0:round(Nvis*PRI*fs)-1)/fs;
A = sqrt(2*db2lin(Pt_dBW - [0 Lt])*k.Z0);      % Tx, echo 1, echo 2 amplitudes
fig = figure('Name', '1.3 Range ambiguity', 'Position', [100 100 900 700]);
labels = {'Transmitted pulses (one per PRI)', ...
          sprintf('Echo target 1 (%d m): %.2f \\mus < PRI, same PRI', R(1), td(1)*1e6), ...
          sprintf('Echo target 2 (%d m): %.2f \\mus = 2 PRI + %.2f \\mus', R(2), td(2)*1e6, mod(td(2),PRI)*1e6)};
delays = [0 td]; cols = {[0.5 0.5 0.5], [0 0.35 0.75], [0.8 0.1 0.1]};
for row = 1:3
    subplot(4,1,row); hold on; grid on;
    e = zeros(size(tv));
    for p = 1:Nvis
        t0 = (p-1)*PRI + delays(row);
        if t0 + tau <= tv(end)
            n0 = ceil(t0*fs - 1e-6);
            e(n0 + (1:round(tau*fs))) = 1;
            text(t0*1e6 + 0.1, 1.25*A(row), sprintf('#%d', p), 'Color', cols{row}, 'FontWeight', 'bold');
        end
        xline((p-1)*PRI*1e6, 'k:');
    end
    area(tv*1e6, A(row)*e, 'FaceColor', cols{row}, 'FaceAlpha', 0.35, 'EdgeColor', cols{row});
    ylim([0 1.6*A(row)]); xlim([0 Nvis*PRI*1e6]); ylabel('V'); title(labels{row});
end
subplot(4,1,4);
[~, ~, e1] = pulse_train(tau, PRI, fc, Pt_dBW, Lt(1), Nvis, fs, R(1));
[~, ~, e2] = pulse_train(tau, PRI, fc, Pt_dBW, Lt(2), Nvis, fs, R(2));
area(tv*1e6, A(2)*e1 + A(3)*e2, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'k'); grid on;
for p = 1:Nvis, xline((p-1)*PRI*1e6, 'k:'); end
xlim([0 Nvis*PRI*1e6]); xlabel('Time [\mus]'); ylabel('V');
title('Receiver (steady state): two echoes per PRI, the far one folded to 352 m');
save_figure(fig, 's01_range_ambiguity');

% Discussion (1.3)
%  * Target 1 (450 m) lies between the blind range (150 m) and Rmax (750 m):
%    its echo returns within its own PRI and is ranged correctly.
%  * Target 2 (1852 m) echoes 12.35 us later, after two more pulses have been
%    sent; the receiver times it from the latest pulse -> ghost at 352 m.
%  * Remedies: lower the PRF (PRI > 12.35 us, PRF < 81 kHz), staggered /
%    multiple PRFs (the true target stays put, ghosts move; Chinese remainder
%    theorem), or pulse-to-pulse coding to identify which pulse an echo belongs to.

%% 1.4 Spectra: single pulse, burst, modulated burst
env1 = double(t < tau);
[f, E1] = spectrum_fft(env1, fs);
[~, E5] = spectrum_fft(env, fs);
[~, V5] = spectrum_fft(v/m.Vpeak, fs);

% [2026] closed-form check. The DFT of a rectangle of M samples is a discrete
% sinc D_M(x) = |sin(pi x M) / sin(pi x)| (-> M tau/T sinc(f tau) for f << fs),
% and repeating it Np times every PRI multiplies by D_Np(f PRI):
%   |X(f)| = D_M(f/fs) D_Np(f PRI) / N,   M = tau fs,  N = record length.
N = numel(t); M = round(tau*fs);
E1_th = dirichlet_mag(f/fs, M) / N;
E5_th = E1_th .* dirichlet_mag(f*PRI, Np);
V5_th = 0.5*(interp1(f, E5_th, f - fc, 'linear', 0) + interp1(f, E5_th, f + fc, 'linear', 0));
err = [max(abs(E1 - E1_th)), max(abs(E5 - E5_th))] ./ [max(E1), max(E5)];
fprintf('\n=== 1.4 Spectra ===\n');
fprintf('sinc nulls every 1/tau = %.1f MHz | lines every PRF = %.0f kHz | line width 1/(Np PRI) = %.0f kHz\n', ...
        1/tau/1e6, m.PRF/1e3, 1/(Np*PRI)/1e3);
fprintf('[2026] max |FFT - closed form| / peak: single pulse %.1e, burst %.1e\n', err);

fig = figure('Name', '1.4 Spectra', 'Position', [100 100 900 700]);
ttl = {'(a) Single 1 \mus pulse: sinc envelope, nulls at k/\tau', ...
       '(b) 5-pulse burst: lines every PRF = 200 kHz under the same sinc', ...
       '(c) Modulated burst: (b) translated to \pm f_c = \pm5 MHz'};
data = {E1, E5, V5}; th = {E1_th, E5_th, V5_th};
for i = 1:3
    subplot(3,1,i);
    plot(f/1e6, data{i}, 'Color', [0 0.35 0.75]); hold on; grid on;
    plot(f/1e6, th{i}, 'r:', 'LineWidth', 1.2);
    xlim([-10 10]); ylabel('|X(f)|'); title(ttl{i});
    if i == 1, legend('FFT', 'Closed form [2026]'); end
end
xlabel('Frequency [MHz]');
save_figure(fig, 's01_spectra');

% Discussion (1.4)
%  (a) A rectangular pulse of width tau has a sinc spectrum: main lobe +-1/tau,
%      nulls at k/tau. Bandwidth ~ 1/tau, so resolution dR = c tau/2 = c/(2B).
%  (b) Repeating the pulse every PRI samples that sinc in frequency: lines every
%      PRF, 5 lines per lobe (1 MHz / 200 kHz). Lines are 1/(Np PRI) = 40 kHz
%      wide because the burst is finite; integrating more pulses narrows them.
%  (c) Multiplying by cos(2 pi fc t) shifts (b) to +-fc with half the amplitude.

%% Local function
function D = dirichlet_mag(x, n)
% |sin(pi x n) / sin(pi x)|, with its limit n where sin(pi x) = 0 (x integer).
s = sin(pi*x);
D = abs(sin(pi*x*n) ./ s);
D(abs(s) < 1e-9) = n;
end
