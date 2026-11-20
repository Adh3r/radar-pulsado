function tests = test_mti
%TEST_MTI  Delay-line cancellers, frequency response, clutter attenuation (Lab 7).
tests = functiontests(localfunctions);
end

function setupOnce(~)
run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'startup.m'));
end

function testCancellerCoefficientsAndOutput(tc)
[y, h] = mti_canceller([1 4 9 16 25], 1);
verifyEqual(tc, h, [1 -1]); verifyEqual(tc, y, [3 5 7 9]);
[y, h] = mti_canceller([1 4 9 16 25].', 2);
verifyEqual(tc, h, [1 -2 1]); verifyEqual(tc, y, [2; 2; 2]);        % 2nd difference of n^2
X = repmat([5 5 5 5], 3, 1);                                         % constant across pulses
verifyEqual(tc, mti_canceller(X, 1), zeros(3, 3));                   % stationary clutter cancelled
verifySize(tc, mti_canceller(rand(10, 32), 2), [10 30]);
end

function testResponseFormula(tc)
% eqs. 7.2-7.3: |H_SC| = 2|sin(w PRI/2)|, |H_DC| = 4 sin^2(w PRI/2)
f = [0 0.1 0.25 0.5 1]; PRI = 1;
verifyEqual(tc, mti_response(f, PRI, 1), 2*abs(sin(pi*f)), 'AbsTol', 1e-12);
verifyEqual(tc, mti_response(f, PRI, 2), 4*sin(pi*f).^2, 'AbsTol', 1e-12);
verifyEqual(tc, mti_response(1e4, 1e-4, 1), 0, 'AbsTol', 1e-12);    % blind speed: fd = PRF
end

function testResponseMatchesFilterOnTones(tc)
% response of the time-domain filter to a complex tone equals the formula
n = (0:199).'; PRI = 1e-4;
for fd = [300 2400 5000]
    for order = 1:2
        y = mti_canceller(exp(1i*2*pi*fd*n*PRI), order);
        verifyEqual(tc, mean(abs(y)), mti_response(fd, PRI, order), 'RelTol', 1e-10);
    end
end
end

function testNoisePowerGain(tc)
rng(31); x = randn(2e5, 1);
verifyEqual(tc, var(mti_canceller(x, 1)), 2, 'RelTol', 0.02);
verifyEqual(tc, var(mti_canceller(x, 2)), 6, 'RelTol', 0.02);
end

function testClutterAttenuation(tc)
% closed form vs small-spread approximations and vs simulation
s = 1e-3;                                                 % sigma_f * PRI
e = 2*pi^2*s^2;
verifyEqual(tc, mti_clutter_attenuation(s, 1, 1), -10*log10(2*(1 - exp(-e))), 'AbsTol', 1e-9);
verifyEqual(tc, mti_clutter_attenuation(s, 1, 2), -10*log10(12*e^2), 'AbsTol', 0.01);
[CA, I] = mti_clutter_attenuation(s, 1, 2);
verifyEqual(tc, I - CA, 10*log10(6), 'AbsTol', 1e-12);
rng(32); c = clutter_slow_time(20000, 4, 0.02, 1, 1);
sim = 10*log10(mean(abs(c(:)).^2) / mean(abs(reshape(mti_canceller(c, 1), [], 1)).^2));
verifyEqual(tc, sim, mti_clutter_attenuation(0.02, 1, 1), 'AbsTol', 0.2);
end

function testClutterSequenceStatistics(tc)
rng(33); sf = 0.05; c = clutter_slow_time(20000, 6, sf, 1, 3);
verifyEqual(tc, mean(abs(c(:)).^2), 3, 'RelTol', 0.02);              % power
r1 = mean(c(:, 2:end) .* conj(c(:, 1:end-1)), 'all') / mean(abs(c(:)).^2);
verifyEqual(tc, real(r1), exp(-2*pi^2*sf^2), 'AbsTol', 0.01);       % lag-1 correlation
c0 = clutter_slow_time(5, 4, 0, 1, 1);                               % sigma_f = 0: static
verifyEqual(tc, c0(:, 1), c0(:, 4), 'AbsTol', 1e-12);
end
