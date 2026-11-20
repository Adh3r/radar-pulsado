function tests = test_cfar
%TEST_CFAR  Cell-averaging CFAR: scaling factor, Pfa model, detector and Pd (Lab 6).
tests = functiontests(localfunctions);
end

function setupOnce(~)
run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'startup.m'));
end

function testLab62ScalingFactor(tc)
% M = 40, Pfa = 0.01 -> alpha = 40 (0.01^(-1/40) - 1)
verifyEqual(tc, ca_cfar_alpha(0.01, 40), 40*(10^(2/40) - 1), 'RelTol', 1e-12);
verifyEqual(tc, ca_cfar_alpha(0.01, 40), 4.8807, 'AbsTol', 1e-4);
end

function testPfaModelRoundTrip(tc)
for M = [8 16 40]
    for N = [1 4 32]
        a = ca_cfar_alpha([1e-2 1e-4 1e-6], M, N);
        verifyEqual(tc, ca_cfar_pfa(a, M, N), [1e-2 1e-4 1e-6], 'RelTol', 1e-8);
    end
end
verifyEqual(tc, ca_cfar_pfa(3, 10), (1 + 3/10)^(-10), 'RelTol', 1e-12);   % eq. 6.1
end

function testLargeMTendsToFixedThreshold(tc)
verifyEqual(tc, ca_cfar_alpha(1e-6, 1e7), -log(1e-6), 'RelTol', 1e-5);
verifyGreaterThan(tc, ca_cfar_alpha(1e-6, 16), -log(1e-6));                 % CFAR needs more
end

function testDetectorWindowAndThreshold(tc)
p = (1:20).';
[det, T, ne] = ca_cfar(p, 4, 1, 2);                % 2 training + 1 guard per side
verifyTrue(tc, all(isnan(ne([1:3, 18:20]))));       % incomplete windows
verifyEqual(tc, ne(10), mean(p([7 8 12 13])), 'AbsTol', 1e-12);
verifyEqual(tc, T(10), 2*ne(10), 'AbsTol', 1e-12);
verifyFalse(tc, any(det([1:3, 18:20])));
[detR, TR] = ca_cfar(p.', 4, 1, 2);                % row input keeps its orientation
verifySize(tc, TR, [1 20]); verifyEqual(tc, detR, det.');
verifyError(tc, @() ca_cfar(p, 5, 1, 2), 'ca_cfar:train');
end

function testMonteCarloPfaIsNoiseIndependent(tc)
% CFAR property: same Pfa whatever the noise power
rng(21); M = 24; a = ca_cfar_alpha(1e-2, M); n = 4e5;
for scale = [1e-15 1 1e6]
    P = scale * -log(rand(M + 7, n));               % exponential noise, any power
    det = ca_cfar(P, M, 3, a);
    verifyEqual(tc, mean(det(M/2 + 4, :)), 1e-2, 'RelTol', 0.05);
end
end

function testMonteCarloIntegratedPfa(tc)
rng(22); M = 16; N = 8; a = ca_cfar_alpha(1e-3, M, N); n = 1e6;
X = reshape(sum(-log(rand(N, (M + 1)*n)), 1), M + 1, n);   % cells = sums of N pulses
cut = X(M/2 + 1, :); tr = mean(X([1:M/2, M/2+2:end], :), 1);
verifyEqual(tc, mean(cut > a*tr), 1e-3, 'RelTol', 0.1);
end

function testSwerling1Pd(tc)
% fixed threshold: Pd = Pfa^(1/(1+SNR)); CA-CFAR reduces to it as M -> inf
verifyEqual(tc, pd_swerling1(10, 1e-6), 1e-6^(1/11), 'RelTol', 1e-12);
verifyEqual(tc, pd_swerling1(10, 1e-6, 1e6), pd_swerling1(10, 1e-6), 'RelTol', 1e-4);
verifyLessThan(tc, pd_swerling1(15, 1e-6, 16), pd_swerling1(15, 1e-6));     % CFAR loss
rng(23); n = 2e5; M = 40; a = ca_cfar_alpha(1e-4, M); snr = db2lin(15);
X = -log(rand(M + 1, n)); X(M/2 + 1, :) = X(M/2 + 1, :) * (1 + snr);
pd = mean(X(M/2 + 1, :) > a * mean(X([1:M/2, M/2+2:end], :), 1));
verifyEqual(tc, pd, pd_swerling1(15, 1e-4, M), 'AbsTol', 0.005);
end
