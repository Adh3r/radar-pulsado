function tests = test_scene
%TEST_SCENE  End-to-end demo scene: reproducibility, calibration and detections.
tests = functiontests(localfunctions);
end

function setupOnce(~)
run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'startup.m'));
end

function testReproducible(tc)
sc = scene_default();
a = simulate_scene(sc); b = simulate_scene(sc);
verifyEqual(tc, a.rx, b.rx);
end

function testCalibration(tc)
sc = scene_default(); out = simulate_scene(sc);
verifySize(tc, out.rx, [round(sc.PRI*sc.fs), sc.n_pulses]);
verifyEqual(tc, mean(abs(out.noise(:)).^2), 1, 'RelTol', 0.02);          % unit noise power
verifyTrue(tc, all(out.rx(out.blind, :) == 0, 'all'));                   % blanking
i = 2;                                                                   % echo power = radar-equation SNR
idx = out.bin(i) + (0:numel(sc.code)-1);
verifyEqual(tc, mean(abs(out.targets(idx, :)).^2, 'all'), db2lin(out.snr_chip_dB(i)), 'RelTol', 1e-10);
verifyEqual(tc, out.range(2) - out.range(1), 15, 'AbsTol', 1e-9);
end

function testTargetsDetectedOutsideClutter(tc)
sc = scene_default(); out = simulate_scene(sc);
K = numel(sc.code);
y = matched_filter(out.rx, sc.code, 'aligned') / sqrt(K);
y(out.blind, :) = 0;
dN = detect_peaks(sum(abs(y).^2, 2), threshold_from_pfa(1e-6, 1, sc.n_pulses));
for i = 2:4                                                              % targets beyond the clutter belt
    verifyTrue(tc, any(abs(dN - out.bin(i)) <= 1), sprintf('target %d not detected', i));
end
% the weak target (12.3 km, 0.1 m^2) is a coin flip for a single pulse
% (theoretical Pd ~ 0.5 after compression) but certain after integration
pd1 = pd_swerling0(out.snr_chip_dB(4) + 10*log10(K), 1e-6);
verifyGreaterThan(tc, pd1, 0.3); verifyLessThan(tc, pd1, 0.7);
end

function testStaticClutterOption(tc)
% sigma_v = 0 gives clutter identical in every pulse; sigma_v > 0 decorrelates it
sc = scene_default(); sc.clutter.sigma_v = 0; out = simulate_scene(sc);
verifyEqual(tc, out.clutter(:, 1), out.clutter(:, end));
sc.clutter.sigma_v = 1; out = simulate_scene(sc);
verifyGreaterThan(tc, max(abs(out.clutter(:, 1) - out.clutter(:, end))), 0);
end

function testMtiRemovesClutterInScene(tc)
% single canceller: residual clutter power close to the closed-form attenuation
sc = scene_default(); sc.targets = zeros(0, 3); out = simulate_scene(sc);
in = out.clutter(out.range > 500 & out.range < 3000, :);
res = mti_canceller(in, 1);
CA = 10*log10(mean(abs(in(:)).^2) / mean(abs(res(:)).^2));
verifyEqual(tc, CA, mti_clutter_attenuation(out.sigma_f, sc.PRI, 1), 'AbsTol', 1);
end
