%% S04 - MATCHED FILTER: Barker pulse compression, code mismatch and Doppler   (Lab 4)
%  Exercises 4.2-4.7 using src/waveform and src/processing, with the metrics
%  the lab asks for (main-to-sidelobe level, number of lobes, output length).
%  Extensions 2026: measured compression gain, Doppler-loss curve and the
%  range-Doppler ambiguity function of Barker-13.

clear; close all;
run(fullfile(fileparts(mfilename('fullpath')), '..', 'startup.m'));
rng(4);
Kset  = [5 7 11 13];
ph    = @(K) 180*(barker_code(K) < 0);          % Barker code -> chip phases [deg]

%% 4.2 Phases of Barker-coded pulses (1 us pulse, 20 samples/chip, SNR 30 dB)
tau = 1e-6; ns = 20; SNR = 30;
fig = figure('Name', '4.2 Barker phases', 'Position', [100 100 900 600]);
for i = 1:4
    K = Kset(i);
    [s, ~, sn, Ts] = phase_coded_pulse(tau, ns, ph(K), [], SNR);
    t = (0:numel(s)-1)*Ts*1e6;
    subplot(2,2,i);
    plot(t, rad2deg(angle(sn)), '.', 'Color', [0.85 0.3 0.3], 'MarkerSize', 4); hold on;
    stairs(t, rad2deg(angle(s)), 'b', 'LineWidth', 1.4); grid on; ylim([-200 200]);
    xlabel('Time [\mus]'); ylabel('Phase [deg]');
    title(sprintf('Barker-%d: chip = %.0f ns', K, tau/K*1e9));
end
legend('With noise (30 dB)', 'Clean', 'Location', 'southwest');
save_figure(fig, 's04_barker_phases');
% Discussion: the pulse length is fixed (1 us), so longer codes have shorter
% chips (200 -> 77 ns) and hence wider bandwidth B ~ K/tau. The phase takes
% only 0 / 180 deg; at 30 dB the noise adds 1/sqrt(2 SNR) = 1.3 deg rms of phase jitter.
% (angle() reports 180 deg as +-180, the same phase.)

%% 4.3 Matched-filter output and sidelobe metrics
fprintf('\n=== 4.3 Matched filter, Barker codes (noise-free metrics) ===\n');
fprintf('%5s %10s %12s %8s %12s %14s\n', 'code', 'PSL [dB]', 'theory [dB]', 'lobes', 'length [us]', 'mainlobe [us]');
fig = figure('Name', '4.3 Matched filter', 'Position', [100 100 1000 750]);
for i = 1:4
    K = Kset(i);
    [s, ~, sn, Ts] = phase_coded_pulse(tau, ns, ph(K), [], SNR);
    y  = matched_filter(sn, s);
    m  = sidelobe_metrics(matched_filter(s, s));
    lag = ((1:numel(y)) - numel(s))*Ts*1e6;
    fprintf('B-%-3d %10.2f %12.2f %8d %12.2f %14.3f\n', K, m.PSL_dB, 20*log10(1/K), ...
            m.n_lobes, m.length*Ts*1e6, m.mainlobe_width*Ts*1e6);
    subplot(4,2,2*i-1); plot(lag, abs(y)/max(abs(y)), 'b'); grid on;
    ylabel('|s_0|'); title(sprintf('Barker-%d (linear)', K));
    subplot(4,2,2*i);   plot(lag, 20*log10(abs(y)/max(abs(y))), 'r'); grid on; ylim([-40 0]);
    yline(20*log10(1/K), 'k--'); ylabel('dB'); title(sprintf('Barker-%d (dB): PSL = %.1f dB', K, m.PSL_dB));
end
subplot(4,2,7); xlabel('Lag [\mus]'); subplot(4,2,8); xlabel('Lag [\mus]');
save_figure(fig, 's04_matched_filter');
% Discussion: peak = pulse energy; every Barker sidelobe has magnitude 1 chip
% (peak K), so PSL = 20 log10(1/K): -14.0, -16.9, -20.8, -22.3 dB. The output
% has K lobes (main + (K-1)/2 per side, odd lags are zero), lasts 2*tau
% (2 us = 2*K*ns - 1 samples) and its main lobe is 2 chips wide null to null:
% the compressed resolution is c*tau_chip/2 instead of c*tau/2, a factor K.

%% 4.4 Barker-13 for SNR from 30 to -10 dB
SNRs = 30:-10:-10;
fig = figure('Name', '4.4 SNR sweep', 'Position', [100 100 900 900]);
fprintf('\n=== 4.4 Barker-13, SNR sweep ===\n');
for i = 1:numel(SNRs)
    [s, ~, sn, Ts] = phase_coded_pulse(tau, ns, ph(13), [], SNRs(i));
    y = matched_filter(sn, s);
    side = abs(y); side(numel(s) + (-2*ns:2*ns)) = 0;
    fprintf('SNR %4d dB: peak/max sidelobe = %5.2f (noise-free 13.0)\n', SNRs(i), max(abs(y))/max(side));
    subplot(numel(SNRs), 2, 2*i-1); plot((0:numel(sn)-1)*Ts*1e6, abs(sn), 'k'); grid on;
    ylabel('|s_n|'); title(sprintf('Input, SNR = %d dB', SNRs(i)));
    subplot(numel(SNRs), 2, 2*i);   plot(((1:numel(y)) - numel(s))*Ts*1e6, abs(y), 'b'); grid on;
    title(sprintf('Matched-filter output, SNR = %d dB', SNRs(i)));
end
subplot(numel(SNRs),2,2*numel(SNRs)-1); xlabel('Time [\mus]');
subplot(numel(SNRs),2,2*numel(SNRs));   xlabel('Lag [\mus]');
save_figure(fig, 's04_snr_sweep');
% Discussion: the peak stays at the pulse energy while the noise floor rises.
% Down to 0 dB the pulse is invisible at the input but the peak is clear at
% the output: the filter integrates 260 samples coherently (24 dB of gain,
% measured below). At -10 dB the output SNR is still ~14 dB, but sidelobes are
% now set by noise rather than by the code; the main-lobe width does not change.

%% [2026] Measured compression gain
n_trials = 4000; SNRin = -10;
[s, ~, ~, Ts] = phase_coded_pulse(tau, ns, ph(13), [], SNRin);
N = numel(s); Pn = 1/db2lin(SNRin);
x = s(:) + sqrt(Pn/2)*(randn(N, n_trials) + 1i*randn(N, n_trials));
y = matched_filter(x, s, 'aligned');           % peak at index 1 (echo starts at sample 1)
yn = matched_filter(sqrt(Pn/2)*(randn(N, n_trials) + 1i*randn(N, n_trials)), s, 'aligned');
gain = lin2db(abs(mean(y(1,:)))^2 / var(yn(1,:))) - SNRin;
fprintf('\n=== [2026] Compression gain, Barker-13 x %d samples/chip ===\n', ns);
fprintf('measured %.2f dB, theory 10 log10(%d) = %.2f dB\n', gain, N, 10*log10(N));

%% 4.5 Barker-5 filter fed with B-5, B-7, B-11 and B-13 (chip 0.1 us)
tchip = 0.1e-6; ns = 20;
s5 = phase_coded_pulse(tchip*5, ns, ph(5), [], SNR);
fig = figure('Name', '4.5 Code mismatch', 'Position', [100 100 1000 420]);
fprintf('\n=== 4.5 Barker-5 filter, mismatched inputs ===\n');
for i = 1:4
    K = Kset(i);
    [~, ~, sn] = phase_coded_pulse(tchip*K, ns, ph(K), [], SNR);
    y = matched_filter(sn, s5);
    subplot(1,2,1); plot(real(y), 'LineWidth', 1.2); hold on;
    subplot(1,2,2); plot(abs(y),  'LineWidth', 1.2); hold on;
    fprintf('input B-%-2d: max |y| = %6.1f  (matched peak = %d)\n', K, max(abs(y)), 5*ns);
end
subplot(1,2,1); grid on; title('Real part'); xlabel('Samples'); ylabel('Re\{s_0\}');
subplot(1,2,2); grid on; title('Magnitude'); xlabel('Samples'); ylabel('|s_0|');
legend(compose('input B-%d', Kset), 'Location', 'northeast');
save_figure(fig, 's04_code_mismatch');
% Discussion: only the B-5 input yields the clean triangle with peak 5 chips
% and 1-chip sidelobes. Longer codes are not replicas of the filter template:
% their outputs are long and irregular. B-11 even contains the negative of
% B-5 (chips 4-8), so it produces a full-height peak (|y| = 100) at the wrong
% lag: a false target. B-7 and B-13 peak at 60 (3 chips).
% The matched filter is a template detector: it is optimal only for the
% waveform it was designed for, which is also what allows code-division use.

%% 4.6 Doppler: Barker-13, f0 = 2.8 GHz, 600 and 60 000 km/h (10 samples/chip)
f0 = 2.8e9; ns = 10;
[s, ~, sn, Ts] = phase_coded_pulse(tchip*13, ns, ph(13), [], SNR);
v_kmh = [0 600 60000];
fig = figure('Name', '4.6 Doppler');
fprintf('\n=== 4.6 Doppler, Barker-13 ===\n');
for v = v_kmh
    fd = doppler_frequency(v/3.6, f0);
    y = matched_filter(apply_doppler(sn, fd, Ts), s);
    plot(((1:numel(y)) - numel(s))*Ts*1e6, abs(y), 'LineWidth', 1.4); hold on;
    fprintf('v = %6d km/h: fd = %7.2f kHz, phase drift over the pulse = %6.1f deg, peak loss = %5.2f dB\n', ...
            v, fd/1e3, 360*fd*numel(s)*Ts, 20*log10(max(abs(y))/(13*ns)));
end
grid on; xlabel('Lag [\mus]'); ylabel('|s_0|');
legend(compose('v_r = %d km/h', v_kmh)); title('Barker-13 matched filter with Doppler shift');
save_figure(fig, 's04_doppler');
% Discussion: at 600 km/h (fd = 3.1 kHz) the phase drifts only 1.5 deg over
% the 1.3 us pulse: the output is unchanged. At 60 000 km/h (fd = 311 kHz)
% it drifts ~146 deg, chips stop adding in phase, the peak drops by 2.5 dB
% and sidelobes rise. Barker codes are Doppler-sensitive: the tolerable fd is
% a fraction of 1/tau (here 770 kHz); beyond it a Doppler filter bank is needed.

%% [2026] Doppler loss curve and ambiguity function
fd_axis = linspace(-2/(13*tchip), 2/(13*tchip), 201);
s1 = repelem(barker_code(13), 4); Ts1 = tchip/4;          % 4 samples/chip
[A, delay] = ambiguity_function(s1, Ts1, fd_axis);
loss_sim = 20*log10(A(:, numel(s1)));                     % zero-lag cut
loss_th  = 20*log10(abs(sin(pi*fd_axis*numel(s1)*Ts1) ./ (numel(s1)*sin(pi*fd_axis*Ts1))));
loss_th(fd_axis == 0) = 0;
fig = figure('Name', '2026 Ambiguity function', 'Position', [100 100 1100 420]);
subplot(1,2,1);
imagesc(delay*1e6, fd_axis/1e3, 20*log10(A + eps)); axis xy; caxis([-30 0]); colorbar;
xlabel('Delay [\mus]'); ylabel('Doppler [kHz]'); title('Barker-13 ambiguity |\chi(\tau, f_d)| [dB]');
subplot(1,2,2);
plot(fd_axis/1e3, loss_sim, 'b', 'LineWidth', 1.4); hold on; grid on;
plot(fd_axis/1e3, loss_th, 'r:', 'LineWidth', 1.6);
xline(doppler_frequency(60000/3.6, f0)/1e3, 'k--', '60 000 km/h');
ylim([-40 1]); xlabel('Doppler [kHz]'); ylabel('Peak loss at zero lag [dB]');
title('Doppler loss'); legend('Simulated', 'Closed form |sin(\pi f_d\tau)/(N sin(\pi f_d T_s))|', 'Location', 'south');
save_figure(fig, 's04_ambiguity');
fprintf('\n=== [2026] Doppler loss: max |sim - closed form| = %.2e dB\n', ...
        max(abs(loss_sim(:).' - loss_th) .* (loss_th > -40)));

%% 4.7 Is there a Barker code of length 6? (chip 0.1 us, 1 sample/chip)
cands = [ones(32,1), 1 - 2*(dec2bin(0:31, 5) - '0')];        % all 2^5 sign-normalised codes
fig = figure('Name', '4.7 Length-6 search'); hold on;
best = inf;
for i = 1:size(cands, 1)
    r = abs(conv(cands(i,:), flip(cands(i,:))));
    plot(-5:5, r, 'Color', [0.75 0.75 0.75]);
    rs = r; rs(6) = 0;
    if max(rs) < best, best = max(rs); ibest = i; end
end
r = abs(conv(cands(ibest,:), flip(cands(ibest,:))));
plot(-5:5, r, 'r-o', 'LineWidth', 1.6); yline(1, 'k--', 'Barker limit');
grid on; xlabel('Lag [chips]'); ylabel('|autocorrelation|');
title(sprintf('All 32 length-6 binary codes: best max sidelobe = %d > 1', best));
save_figure(fig, 's04_barker6_search');
fprintf('\n=== 4.7 Exhaustive search ===\n');
for K = [5 6 7 13]
    fprintf('K = %2d: %d code(s) with Barker sidelobes\n', K, size(barker_search(K), 1));
end
% Discussion: none of the 32 sign-normalised length-6 codes (64 counting
% negatives) keeps every sidelobe <= 1: the best reaches 2. The same search
% finds the known codes for K = 5, 7 and 13 (positive controls; each found
% twice, the code and its time reversal), so the negative result for K = 6 is
% not an artefact of the method.

%% Exploration: signal anatomy (from the study demos written for Lab 4)
% Why take |s0| at all, and what does Doppler do to the samples? These views
% come from the exploratory demos written while studying the lab: the BPSK
% waveform on a (low, visual) carrier, its baseband I component, the I/Q
% cloud, and the phase ramp / constellation arc a Doppler shift imposes.
rng(5);
f_vis = 10e6;                                            % visual carrier only
[s, ~, sn, Ts] = phase_coded_pulse(1e-6, 50, ph(5), [], 15);
t = (0:numel(s)-1)*Ts;
fd = doppler_frequency(60000/3.6, 2.8e9);
s13 = phase_coded_pulse(13*tchip, 20, ph(13), [], 30);
[~, ~, sn13, Ts13] = phase_coded_pulse(13*tchip, 20, ph(13), [], 30);
sd13 = apply_doppler(sn13, fd, Ts13);
t13 = (0:numel(sn13)-1)*Ts13;
fig = figure('Name', 'Exploration: signal anatomy', 'Position', [100 100 1000 700]);
subplot(2,2,1);
plot(t*1e6, real(s.*exp(1i*2*pi*f_vis*t)), 'b'); hold on;
plot(t*1e6, real(s), 'r', 'LineWidth', 1.6); grid on; ylim([-1.5 1.5]);
xlabel('Time [\mus]'); title('Barker-5 BPSK: a -1 chip inverts the carrier');
legend('RF (10 MHz, visual)', 'baseband I', 'Location', 'southwest');
subplot(2,2,2);
plot(real(sn), imag(sn), '.', 'Color', [0.85 0.3 0.3], 'MarkerSize', 5); hold on;
plot([1 -1], [0 0], 'bx', 'MarkerSize', 12, 'LineWidth', 2.5); grid on; axis equal;
xlim([-2 2]); ylim([-2 2]); xlabel('I'); ylabel('Q');
title('I/Q at 15 dB: the two BPSK points become noise clouds');
subplot(2,2,3);
plot(t13*1e6, mod(rad2deg(angle(sn13)) + 90, 360) - 90, 'b.', 'MarkerSize', 4); hold on;
plot(t13*1e6, mod(rad2deg(angle(sd13)) + 90, 360) - 90, 'r.', 'MarkerSize', 4); grid on;
xlabel('Time [\mus]'); ylabel('Phase [deg]'); ylim([-90 270]);
title(sprintf('Barker-13 phase: static vs 60 000 km/h (f_d = %.0f kHz ramp)', fd/1e3));
legend('static', 'Doppler', 'Location', 'southwest');
subplot(2,2,4);
plot(real(sn13), imag(sn13), 'b.', real(sd13), imag(sd13), 'r.', 'MarkerSize', 6); hold on;
th = linspace(0, 2*pi*fd*t13(end), 100); plot(cos(th), sin(th), 'k:', 'LineWidth', 1.3);
grid on; axis equal; xlim([-1.5 1.5]); ylim([-1.5 1.5]); xlabel('I'); ylabel('Q');
title('Doppler rotates the constellation along an arc');
save_figure(fig, 's04_signal_anatomy');
% The matched-filter output is complex: its phase depends on the unknown echo
% phase (and on Doppler), so detection uses |s0| (envelope), which is
% insensitive to that phase. The arc shows why chips stop adding coherently
% when the Doppler phase drift over the pulse becomes a sizeable fraction of 2*pi.

%% Exploration: Barker-13 filter fed with shorter codes (fixed 13-chip window)
% The converse of 4.5: a filter matched to B-13 receives B-5, B-7, B-11 and
% B-13, each padded with noise to the same 13-chip window.
ns = 20; Pn = 1/db2lin(SNR);
fig = figure('Name', 'Exploration: Barker-13 filter', 'Position', [100 100 1000 420]);
fprintf('\n=== Exploration: Barker-13 filter, shorter inputs ===\n');
for i = 1:4
    K = Kset(i);
    [~, ~, sn] = phase_coded_pulse(tchip*K, ns, ph(K), [], SNR);
    tail = sqrt(Pn/2)*(randn(1, (13-K)*ns) + 1i*randn(1, (13-K)*ns));
    y = matched_filter([sn, tail], s13);
    subplot(1,2,1); plot(real(y), 'LineWidth', 1.2); hold on;
    subplot(1,2,2); plot(abs(y),  'LineWidth', 1.2); hold on;
    fprintf('input B-%-2d: max |y| = %6.1f (matched peak = %d)\n', K, max(abs(y)), 13*ns);
end
subplot(1,2,1); grid on; title('Real part (sign of the correlation)'); xlabel('Samples');
subplot(1,2,2); grid on; title('Magnitude |s_0| (what the detector sees)'); xlabel('Samples');
legend(compose('input B-%d', Kset), 'Location', 'northeast');
save_figure(fig, 's04_b13_filter_mismatch');
% The real part keeps the sign: a code segment that matches the INVERTED
% template gives a large negative correlation, which |s0| turns into an
% equally large positive peak. That is how B-11 (which contains -B-5)
% produces a full-height false peak in 4.5: the envelope detector cannot
% tell a match from an inverted match.
