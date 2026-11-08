function tests = test_waveform
%TEST_WAVEFORM  Pulse metrics, pulse train and spectrum (Lab 1).
tests = functiontests(localfunctions);
end

function setupOnce(~)
run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'startup.m'));
end

function testLab12Metrics(tc)
% Lab 1.2 values (identical to the original lab code)
m = pulse_metrics(1e-6, 5e-6, 60, 70);
verifyEqual(tc, [m.PRF m.Rmin m.Rmax m.dR m.duty], [200e3 150 750 150 20], 'AbsTol', 1e-9);
verifyEqual(tc, m.Pavg_dBm, 20 - 10*log10(5), 'AbsTol', 1e-12);   % peak +20 dBm x 20 % duty
verifyEqual(tc, m.Epulse, 0.1*1e-6, 'RelTol', 1e-12);
verifyEqual(tc, m.Vpeak, sqrt(2*0.1*50), 'RelTol', 1e-12);
end

function testEveryPulseHasExactlyMSamples(tc)
% regression for the floating-point boundary bug of the original code
[~, ~, env] = pulse_train(1e-6, 5e-6, 5e6, 60, 70, 5, 100e6, 0);
d = diff([0 env 0]);
verifyEqual(tc, find(d == -1) - find(d == 1), 100*ones(1, 5));
verifyEqual(tc, find(d == 1), 1 + (0:4)*500);
end

function testRangeAmbiguityFolding(tc)
% 1852 m with Rmax = 750 m: echo folds to 2.35 us (352 m)
fs = 100e6;
[~, ~, env] = pulse_train(1e-6, 5e-6, 5e6, 60, 6, 5, fs, 1852);
first = find(env, 1) - 1;                           % 0-based sample of first echo
verifyEqual(tc, first/fs, mod(2*1852/3e8, 5e-6), 'AbsTol', 1/fs);
end

function testSpectrumMatchesClosedForm(tc)
fs = 100e6; [~, v, t, env] = pulsed_radar(1e-6, 5e-6, 5e6, 60, 70, 5, fs, 0); %#ok<ASGLU>
[f, X] = spectrum_fft(env, fs);
N = numel(env);
D = @(x, n) local_dirichlet(x, n);
verifyLessThan(tc, max(abs(X - D(f/fs, 100)/N .* D(f*5e-6, 5))), 1e-12);
verifyEqual(tc, sum(X.^2), sum(env.^2)/N^2*N, 'RelTol', 1e-12);   % Parseval
end

function D = local_dirichlet(x, n)
s = sin(pi*x); D = abs(sin(pi*x*n)./s); D(abs(s) < 1e-9) = n;
end
