function tests = test_processing
%TEST_PROCESSING  Barker codes, coded pulse, matched filter, sidelobes, Doppler, ambiguity (Lab 4).
tests = functiontests(localfunctions);
end

function setupOnce(~)
run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'startup.m'));
end

function testBarkerSidelobeMetrics(tc)
% Lab 4.3: PSL = 20 log10(1/K), K lobes, length 2N-1, main lobe 2 chips
ns = 20;
for K = [5 7 11 13]
    s = repelem(barker_code(K), ns);
    m = sidelobe_metrics(matched_filter(s, s));
    verifyEqual(tc, m.PSL_dB, 20*log10(1/K), 'AbsTol', 1e-9);
    verifyEqual(tc, m.n_lobes, K);
    verifyEqual(tc, m.length, 2*K*ns - 1);
    verifyEqual(tc, m.mainlobe_width, 2*ns);
    verifyEqual(tc, m.peak, K*ns, 'AbsTol', 1e-9);
end
end

function testMatchedFilterModes(tc)
code = barker_code(13).';
x = zeros(100, 3); x(41:53, :) = repmat(code, 1, 3);
y = matched_filter(x, code, 'aligned');
verifySize(tc, y, [100 3]);
[~, ip] = max(abs(y));
verifyEqual(tc, ip, [41 41 41]);                    % peak at the echo start
yf = matched_filter(x(:,1).', code);                % row vector, full
verifySize(tc, yf, [1 112]);
[~, ipf] = max(abs(yf)); verifyEqual(tc, ipf, 41 + 13 - 1);
end

function testBarkerSearch(tc)
verifyEmpty(tc, barker_search(6));                  % Lab 4.7: no length-6 Barker code
for K = [5 7 13]                                     % positive controls
    codes = barker_search(K);
    verifyEqual(tc, size(codes, 1), 2);              % code and its time reversal
    verifyTrue(tc, ismember(barker_code(K), codes, 'rows') || ...
                   ismember(-barker_code(K), codes, 'rows'));
end
end

function testDopplerShiftAndLoss(tc)
verifyEqual(tc, doppler_frequency(60000/3.6, 2.8e9), 2*60000/3.6*2.8e9/3e8, 'RelTol', 1e-12);
s = repelem(barker_code(13), 10); Ts = 1e-8;
fd = 311.11e3;
y = matched_filter(apply_doppler(s, fd, Ts), s);
N = numel(s);
loss_th = abs(sin(pi*fd*N*Ts)/(N*sin(pi*fd*Ts)));
verifyEqual(tc, abs(y(N))/N, loss_th, 'RelTol', 1e-10);   % zero-lag peak loss
end

function testAmbiguityFunction(tc)
s = repelem(barker_code(13), 2); Ts = 0.05e-6;
[A, delay] = ambiguity_function(s, Ts, [0 1e5]);
verifyEqual(tc, max(A(1,:)), 1, 'AbsTol', 1e-12);
verifyEqual(tc, delay(numel(s)), 0, 'AbsTol', 1e-15);
verifyLessThan(tc, A(2, numel(s)), 1);
end

function testCompressionGain(tc)
rng(11); code = barker_code(13); K = 13; n = 2e4;
noise = sqrt(1/2)*(randn(K, n) + 1i*randn(K, n));
y = matched_filter(noise, code, 'aligned');
verifyEqual(tc, var(y(1,:)), K, 'RelTol', 0.03);    % output noise power = energy * sigma^2
end

function testDetectPeaks(tc)
p = [0 5 6 5 0 0 9 0 2];
[idx, val] = detect_peaks(p, 4);
verifyEqual(tc, idx.', [3 7]);
verifyEqual(tc, val.', [6 9]);
end

function testBarkerAutocorrelation(tc)
for K = [2 3 4 5 7 11 13]
    c = barker_code(K);
    r = conv(c, flip(c));
    verifyEqual(tc, r(K), K);
    r(K) = 0;
    verifyLessThanOrEqual(tc, max(abs(r)), 1);
end
verifyError(tc, @() barker_code(6), 'barker_code:length');
end

function testPhaseCodedPulseSnrAndSampling(tc)
rng(1);
[s, n, sn, Ts] = phase_coded_pulse(1.3e-6, 20, 180*(barker_code(13) < 0), [], 10);
verifyEqual(tc, Ts, 0.1e-6/20, 'RelTol', 1e-12);
verifyEqual(tc, numel(s), 260);
verifyEqual(tc, real(s(1:20:end)), barker_code(13), 'AbsTol', 1e-12);
verifyEqual(tc, sn, s + n);
rng(2);                                             % SNR over many samples
[s, n] = phase_coded_pulse(1e-3, 100, zeros(1, 1000), [], 10);
verifyEqual(tc, 10*log10(mean(abs(s).^2)/mean(abs(n).^2)), 10, 'AbsTol', 0.05);
end
