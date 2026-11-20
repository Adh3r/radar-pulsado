function out = simulate_scene(sc)
%SIMULATE_SCENE  Baseband data cube (fast time x slow time) for a scenario.  [2026]
%   out = simulate_scene(sc)     sc from SCENE_DEFAULT
%   Fields of out:
%     rx          complex samples, n_fast x n_pulses (noise power = 1 per sample)
%     noise, clutter, targets   the three components of rx (same size)
%     range       range of each fast-time sample [m]
%     snr_chip_dB per-target single-sample SNR from the radar equation
%     fd          per-target Doppler shift [Hz]
%     bin         per-target fast-time index where its echo starts
%     blind       logical mask of samples blanked while transmitting
%
%   Model
%   * Noise: complex white Gaussian, unit power (radar-equation SNRs are
%     therefore directly the echo powers).
%   * Targets: constant amplitude sqrt(SNR), random initial phase, Doppler
%     phase progression exp(j 2 pi fd (t_fast + m PRI)) across the CPI.
%   * Clutter: one complex Gaussian scatterer per range bin (Rayleigh
%     amplitude) between r_min and r_max, power decaying as R^-exponent,
%     centred at zero Doppler: identical in every pulse if sc.clutter.sigma_v
%     is 0, otherwise decorrelating from pulse to pulse with a Gaussian
%     Doppler spectrum (eq. 7.1 of Lab 7). Each scatterer is illuminated
%     by the full coded pulse, so clutter is the code convolved with the
%     reflectivity profile.
%   * Samples within the transmitted pulse length are blanked (blind range).
%   The random generator is seeded with sc.seed for reproducibility.
rng(sc.seed);
k  = radar_constants();
K  = numel(sc.code);
nf = round(sc.PRI*sc.fs);
np = sc.n_pulses;
t_fast = (0:nf-1).'/sc.fs;
out.range = k.c*t_fast/2;
B = sc.fs;                                         % noise bandwidth ~ 1/tchip

% --- noise ----------------------------------------------------------------
out.noise = sqrt(1/2)*(randn(nf, np) + 1i*randn(nf, np));

% --- clutter (static) -----------------------------------------------------
c = sc.clutter;
inC = out.range >= c.r_min & out.range <= c.r_max;
p_c = zeros(nf, 1);
p_c(inC) = db2lin(c.CNR0_dB) * (c.r_min ./ out.range(inC)).^c.exponent / K;   % per-scatterer power
sigma_v = 0; if isfield(c, 'sigma_v'), sigma_v = c.sigma_v; end
if sigma_v == 0                                    % perfectly static clutter
    refl = sqrt(p_c/2) .* (randn(nf,1) + 1i*randn(nf,1));
    cl = conv(refl, sc.code(:));
    out.clutter = repmat(cl(1:nf), 1, np);
else                                               % [2026] internal motion: Gaussian Doppler spread
    out.sigma_f = 2*sigma_v*sc.fc/k.c;             % velocity spread -> Doppler spread [Hz]
    refl = clutter_slow_time(nf, np, out.sigma_f, sc.PRI, p_c);
    cl = conv2(refl, sc.code(:));                  % each pulse illuminates the scatterers
    out.clutter = cl(1:nf, :);
end

% --- targets --------------------------------------------------------------
nt = size(sc.targets, 1);
out.targets = zeros(nf, np);
out.snr_chip_dB = radar_snr(sc.Pt, sc.G_dBi, sc.fc, B, sc.NF_dB, sc.L_dB, ...
                            sc.targets(:,2), sc.targets(:,1));
out.fd  = doppler_frequency(sc.targets(:,3), sc.fc);
out.bin = round(2*sc.targets(:,1)/k.c*sc.fs) + 1;
slow = (0:np-1)*sc.PRI;
for i = 1:nt
    n0 = out.bin(i);
    if n0 + K - 1 > nf, continue, end              % beyond the unambiguous window
    a  = sqrt(db2lin(out.snr_chip_dB(i))) * exp(1i*2*pi*rand);
    idx = n0:n0+K-1;
    pulse = sc.code(:) .* exp(1i*2*pi*out.fd(i)*t_fast(idx));
    out.targets(idx, :) = out.targets(idx, :) + a * pulse .* exp(1i*2*pi*out.fd(i)*slow);
end

% --- receiver blanking during transmission -------------------------------
out.blind = (1:nf).' <= K;
out.rx = out.noise + out.clutter + out.targets;
out.rx(out.blind, :) = 0;
end
