function tests = test_detection
%TEST_DETECTION  False-alarm time, Albersheim, Marcum Q, thresholds and ROC.
tests = functiontests(localfunctions);
end

function setupOnce(~)
run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'startup.m'));
end

function testFalseAlarmTime(tc)
[Tfa, Pfa] = false_alarm_time([0 10], [1 1e6]);
verifyEqual(tc, Pfa, exp(-[1; 10]), 'RelTol', 1e-12);
verifyEqual(tc, Tfa, exp([1; 10]) ./ [1 1e6], 'RelTol', 1e-12);
end

function testAlbersheimReference(tc)
% Skolnik reference point and the shape of the output
verifyEqual(tc, albersheim_snr(0.9, 1e-6), 13.14, 'AbsTol', 0.01);
verifySize(tc, albersheim_snr([0.5 0.9 0.99], 10.^(-8:-3)), [3 6]);
verifyWarning(tc, @() albersheim_snr(0.1, 1e-6), 'albersheim_snr:range');
end

function testMarcumQKnownValues(tc)
verifyEqual(tc, marcum_q(0, [0.5 2 3 7.4]), exp(-[0.5 2 3 7.4].^2/2), 'RelTol', 1e-8);  % Rayleigh tail
verifyEqual(tc, marcum_q(2, 0), 1, 'AbsTol', 1e-10);
a = 1.5; b = 2.5;                                                     % symmetry identity
verifyEqual(tc, marcum_q(a, b) + marcum_q(b, a), 1 + exp(-(a^2 + b^2)/2)*besseli(0, a*b), 'RelTol', 1e-8);
end

function testAlbersheimAgainstExact(tc)
% Albersheim's empirical fit vs exact Pd over the lab range
Pd = [0.3 0.5 0.9 0.99]; Pfa = 10.^(-10:-3);
S  = albersheim_snr(Pd, Pfa);
err = zeros(size(S));
for i = 1:numel(Pd)
    for j = 1:numel(Pfa)
        err(i,j) = fzero(@(s) pd_swerling0(s, Pfa(j)) - Pd(i), S(i,j)) - S(i,j);
    end
end
verifyLessThan(tc, max(abs(err(:))), 0.3);     % dB
end

function testThresholdsGiveDesignPfa(tc)
verifyEqual(tc, threshold_from_pfa(1e-6, 1), -log(1e-6), 'RelTol', 1e-12);
verifyEqual(tc, threshold_from_pfa(1e-3, 2, 1), -2*log(1e-3), 'RelTol', 1e-12);
% Monte Carlo: n = 1 and n = 8 integrated pulses
rng(7); n = 2e6;
x = (randn(n, 8).^2 + randn(n, 8).^2)/2;      % exponential, unit mean
verifyEqual(tc, mean(x(:,1) > threshold_from_pfa(1e-3, 1)), 1e-3, 'RelTol', 0.1);
verifyEqual(tc, mean(sum(x, 2) > threshold_from_pfa(1e-3, 1, 8)), 1e-3, 'RelTol', 0.1);
end

function testMonteCarloPdMatchesMarcum(tc)
rng(8); n = 2e5; snr = 12; Pfa = 1e-4;
x = sqrt(db2lin(snr))*exp(1i*2*pi*rand(n,1)) + sqrt(1/2)*(randn(n,1) + 1i*randn(n,1));
verifyEqual(tc, mean(abs(x).^2 > threshold_from_pfa(Pfa, 1)), pd_swerling0(snr, Pfa), 'AbsTol', 0.005);
end

function testRocEmpirical(tc)
[pf, pd] = roc_empirical([2 3 4], [0 1 2 3], [-1 1.5 3.5 5]);
verifyEqual(tc, pd,  [1 1 1/3 0], 'AbsTol', 1e-12);
verifyEqual(tc, pf,  [1 0.5 0 0], 'AbsTol', 1e-12);
end
