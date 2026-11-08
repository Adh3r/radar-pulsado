function tests = test_link_budget
%TEST_LINK_BUDGET  Friis cascade, radar equation and detection range.
tests = functiontests(localfunctions);
end

function setupOnce(~)
run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'startup.m'));
end

function testLab2NoiseFigures(tc)
% Values of the original lab code: 1.2086 dB (LNA first), 3.9043 dB (cable first)
Lc = 2.7;
NFA = cascade_noise_figure([37 -Lc -1 -6 -1 11], [1.2 Lc 1 5.5 1 1.1]);
NFB = cascade_noise_figure([-Lc 37 -1 -6 -1 11], [Lc 1.2 1 5.5 1 1.1]);
verifyEqual(tc, [NFA NFB], [1.2086 3.9043], 'AbsTol', 1e-4);
end

function testFriisSpecialCases(tc)
verifyEqual(tc, cascade_noise_figure(20, 3), 3, 'AbsTol', 1e-12);              % single stage
verifyEqual(tc, cascade_noise_figure([-3 -4], [3 4]), 7, 'AbsTol', 1e-12);      % passive chain: losses add
verifyError(tc, @() cascade_noise_figure([1 2], 3), 'cascade_noise_figure:size');
end

function testRadarEquationLaws(tc)
args = {1.5e6, 30, 1.1e9, 50e6, 1.2086, 0};
S100 = radar_snr(args{:}, 1, 100e3);
verifyEqual(tc, S100, 3.278, 'AbsTol', 1e-3);                                   % lab value
verifyEqual(tc, radar_snr(args{:}, 1, 200e3) - S100, -40*log10(2), 'AbsTol', 1e-10);  % R^-4
verifyEqual(tc, radar_snr(args{:}, 0.1, 100e3) - S100, -10, 'AbsTol', 1e-10);    % RCS linear
verifyEqual(tc, radar_snr(args{1:5}, 10, 1, 100e3) - S100, -10, 'AbsTol', 1e-10); % losses
end

function testMaxRangeInvertsRadarEquation(tc)
args = {1.5e6, 30, 1.1e9, 50e6, 1.2086, 0, 1};
R = max_detection_range(args{:}, [5 13.14 20]);
verifyEqual(tc, radar_snr(args{:}, R), [5 13.14 20], 'AbsTol', 1e-10);
end
