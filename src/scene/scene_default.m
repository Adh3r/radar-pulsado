function sc = scene_default()
%SCENE_DEFAULT  Reference scenario for the end-to-end demo.  [2026]
%   sc = scene_default() returns a struct with the radar, waveform,
%   clutter and target definitions. Edit the fields and pass the struct to
%   SIMULATE_SCENE.
%
%   S-band surveillance radar, Barker-13 phase-coded pulse, a coherent
%   processing interval (CPI) of 32 pulses. The receiver noise figure is the
%   LNA-first chain of Lab 2 (scenario A).

% --- radar --------------------------------------------------------------
sc.fc     = 3e9;            % carrier [Hz]
sc.Pt     = 5e3;            % peak power [W]
sc.G_dBi  = 30;             % antenna gain (Tx = Rx)
sc.L_dB   = 3;              % system losses
sc.NF_dB  = cascade_noise_figure([37 -2.7 -1 -6 -1 11], [1.2 2.7 1 5.5 1 1.1]);

% --- waveform -------------------------------------------------------------
sc.code    = barker_code(13);   % chip values +1/-1
sc.tchip   = 0.1e-6;            % chip duration [s] -> 15 m compressed resolution
sc.fs      = 1/sc.tchip;        % 1 sample per chip: fast-time bin = range bin
sc.PRI     = 100e-6;            % -> 15 km unambiguous range
sc.n_pulses = 32;               % CPI length (slow time)

% --- clutter: ground returns, Rayleigh amplitude, ~R^-3 power law ---------
sc.clutter.r_min   = 300;       % [m]
sc.clutter.r_max   = 3500;      % [m]
sc.clutter.CNR0_dB = 45;        % clutter-to-noise ratio per bin at r_min (before compression)
sc.clutter.exponent = 3;
sc.clutter.sigma_v  = 1.0;      % [2026] rms internal-motion velocity [m/s]: illustrative strong-wind value
                                % Gaussian Doppler spread sigma_f = 2 sigma_v / lambda (0 = static)

% --- targets (non-fluctuating, Swerling 0) -------------------------------
%                 range [m]  RCS [m^2]  radial velocity [m/s] (+ approaching)
sc.targets = [      2400      0.01       15     % small drone buried in clutter
                    5500      0.5       120
                    8700      2.0       -60
                   12300      0.1       250 ];  % weak: needs pulse compression
sc.seed = 2026;
end
