%% S02 - RADAR EQUATION: receiver noise figure and SNR vs range   (Lab 2)
%  Exercises 2.1-2.6 using src/link_budget. The 2026 detection-range
%  calculator, which needs Albersheim's equation (Lab 3), is in s03.

clear; close all;
run(fullfile(fileparts(mfilename('fullpath')), '..', 'startup.m'));
k = radar_constants();

%% 2.1 Noise figure of the two receiver chains (Fig. 1 of the lab)
% Passive stage (cable, filter, passive mixer): gain = -loss, NF = +loss.
L_coax = 0.27*10;                                    % RG-8x, 0.27 dB/m x 10 m
stages = {'LNA', 'Coax', 'RF BPF', 'Mixer', 'IF BPF', 'IF amp'};
G_A  = [37, -L_coax, -1, -6, -1, 11];  NF_A = [1.2, L_coax, 1, 5.5, 1, 1.1];   % A: LNA first
G_B  = G_A([2 1 3:end]);               NF_B = NF_A([2 1 3:end]);               % B: cable first
[NFA, ~, cA] = cascade_noise_figure(G_A, NF_A);
[NFB, ~, cB] = cascade_noise_figure(G_B, NF_B);

fprintf('\n=== 2.1 Cascade noise figure ===\n');
fprintf('A (LNA first):   NF = %.3f dB\nB (cable first): NF = %.3f dB  (+%.2f dB)\n', NFA, NFB, NFB - NFA);
fprintf('%-8s %12s %12s\n', 'stage', 'contrib. A', 'contrib. B');
namesB = stages([2 1 3:end]);
for i = 1:numel(stages)
    fprintf('%-8s %12.5f %12.5f   (%s / %s)\n', stages{i}, cA(i), cB(i), stages{i}, namesB{i});
end
% Discussion: in A the LNA's 37 dB divide every later contribution by ~5000,
% so NF ~ NF_LNA. In B the 2.7 dB of cable come before any gain and add
% almost one-to-one to the NF: any loss ahead of the first amplifier adds
% directly to the system NF, which is why LNAs are mounted at the antenna.

%% 2.3 SNR vs range for both chains
Pt = 1.5e6; fc = 1.1e9; B = 50e6; sigma = 1;
G_ant = 30;   % ASSUMPTION: Fig. 1 does not give the antenna gain (37 dB is the
              % LNA). 30 dBi is typical of a 1.1 GHz dish; change it here.
R = (20:0.5:200)*1e3;
SNR_A = radar_snr(Pt, G_ant, fc, B, NFA, 0, sigma, R);
SNR_B = radar_snr(Pt, G_ant, fc, B, NFB, 0, sigma, R);

fprintf('\n=== 2.3 SNR vs range (G = %d dBi) ===\n', G_ant);
fprintf('noise power kT0BF = %.1f dBm (chain A)\n', lin2db(k.kB*k.T0*B*db2lin(NFA)*1e3));
for Rq = [20 50 100 200]*1e3
    fprintf('R = %3.0f km: SNR_A = %6.2f dB, SNR_B = %6.2f dB\n', Rq/1e3, ...
            radar_snr(Pt, G_ant, fc, B, NFA, 0, sigma, Rq), radar_snr(Pt, G_ant, fc, B, NFB, 0, sigma, Rq));
end

fig = figure('Name', '2.3 SNR vs range', 'Position', [100 100 1000 420]);
subplot(1,2,1);
plot(R/1e3, SNR_A, 'b-', R/1e3, SNR_B, 'r--', 'LineWidth', 1.5); grid on; yline(0, 'k:');
xlabel('Range [km]'); ylabel('SNR [dB]'); title('Linear range axis');
legend(sprintf('A: LNA first (NF %.2f dB)', NFA), sprintf('B: cable first (NF %.2f dB)', NFB));
subplot(1,2,2);
semilogx(R/1e3, SNR_A, 'b-', R/1e3, SNR_B, 'r--', 'LineWidth', 1.5); grid on; yline(0, 'k:');
xlabel('Range [km] (log)'); ylabel('SNR [dB]'); title('Log range axis: straight line, -40 dB/decade');
save_figure(fig, 's02_snr_vs_range');
% Discussion: R^-4 becomes a straight line of slope -40 dB/decade (-12 dB per
% range doubling) on semilogx. The NF only shifts the curve: the 2.7 dB gap is
% constant and costs 10^(-2.7/40) -> 14 % of detection range.

%% 2.4-2.6 Sensitivity to RCS, transmitted power and extra losses
cases = struct( ...
    'name',  {'RCS', 'Peak power', 'Extra losses'}, ...
    'vals',  {[0 -10 -20], [0 3 6], [0 10]}, ...
    'leg',   {{'\sigma = 0 dBsm', '\sigma = -10 dBsm', '\sigma = -20 dBsm'}, ...
              {'P_t = 1.5 MW', 'P_t + 3 dB', 'P_t + 6 dB'}, {'L = 0 dB', 'L = 10 dB'}});
fig = figure('Name', '2.4-2.6 Sensitivity', 'Position', [100 100 1200 380]);
fprintf('\n=== 2.4-2.6 Sensitivity (chain A, SNR at 100 km) ===\n');
for c = 1:3
    subplot(1,3,c);
    for v = cases(c).vals
        switch c
            case 1, S = radar_snr(Pt, G_ant, fc, B, NFA, 0, db2lin(v), R);
            case 2, S = radar_snr(Pt*db2lin(v), G_ant, fc, B, NFA, 0, sigma, R);
            case 3, S = radar_snr(Pt, G_ant, fc, B, NFA, v, sigma, R);
        end
        semilogx(R/1e3, S, 'LineWidth', 1.5); hold on;
        fprintf('%-12s %+4d dB -> %6.2f dB  (range factor %.2f)\n', cases(c).name, v, ...
                interp1(R, S, 100e3), 10^((interp1(R, S, 100e3) - interp1(R, SNR_A, 100e3))/40));
    end
    grid on; yline(0, 'k:'); xlabel('Range [km]'); ylabel('SNR [dB]');
    title(cases(c).name); legend(cases(c).leg, 'Location', 'southwest');
end
save_figure(fig, 's02_sensitivity');
% Discussion: Pt, sigma, L (and F) enter the equation linearly, so they only
% shift the curve; only R changes its slope. Because range goes as the 4th
% root: -10 dBsm of RCS -> 0.56x range (the case for stealth), +3/+6 dB of
% power -> only +19/+41 % range (doubling range needs 16x power), 10 dB of
% losses -> 0.56x range. Realistic budgets keep several dB of margin.
