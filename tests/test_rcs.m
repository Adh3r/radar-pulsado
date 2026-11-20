function tests = test_rcs
%TEST_RCS  Sphere (Mie), circular plate (PO), multi-scatterer RCS and Swerling models (Lab 5).
tests = functiontests(localfunctions);
end

function setupOnce(~)
run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'startup.m'));
end

function testSphereLimits(tc)
verifyEqual(tc, rcs_sphere_mie(0.02), 9*0.02^4, 'RelTol', 1e-3);       % Rayleigh: 9 (ka)^4
verifyEqual(tc, mean(rcs_sphere_mie(20:0.01:40)), 1, 'AbsTol', 0.01);   % optical: pi a^2
[pk, i] = max(rcs_sphere_mie(0.5:0.001:2));
verifyEqual(tc, pk, 3.65, 'AbsTol', 0.02);                               % first resonance
verifyEqual(tc, 0.5 + (i-1)*0.001, 1.03, 'AbsTol', 0.02);
end

function testPlateNormalIncidenceAndNull(tc)
r = 0.25; f = 3e9; lambda = 3e8/f;
verifyEqual(tc, rcs_circular_plate(r, f, 0), 4*pi*(pi*r^2)^2/lambda^2, 'RelTol', 1e-12);
verifyEqual(tc, rcs_circular_plate(r, f, 90), 0, 'AbsTol', 1e-20);      % edge-on (PO)
verifyEqual(tc, rcs_circular_plate(r, f, asind(3.8317059702*lambda/(4*pi*r))), 0, 'AbsTol', 1e-8);
end

function testSingleScattererIsIsotropic(tc)
radar = [1e3*ones(37,1), 90*ones(37,1), (0:10:360)'];
verifyEqual(tc, rcs_scatterers(3, [0.4 -0.2 0], radar, 3e9), 3*ones(37,1), 'AbsTol', 1e-9);
end

function testTwoSpheresClosedForm(tc)
fc = 1.5e9; kw = 2*pi*fc/3e8; L = 0.37; phi = (0:1:360)';
sig = db2lin(rcs_scatterers([0; 0], [-L/2 0 0; L/2 0 0], [1e5*ones(size(phi)), 90*ones(size(phi)), phi], fc));
verifyEqual(tc, sig, 4*cos(kw*L*cosd(phi)).^2, 'AbsTol', 1e-3);         % far field
verifyEqual(tc, max(sig), 4, 'AbsTol', 1e-6);                            % coherent sum: +6 dBsm
end

function testSwerlingPdfs(tc)
x = linspace(0, 50, 200001);
for model = [1 3]
    p = swerling_pdf(x, 2, model);
    verifyEqual(tc, trapz(x, p), 1, 'AbsTol', 1e-6);                     % normalised
    verifyEqual(tc, trapz(x, x.*p), 2, 'AbsTol', 1e-5);                  % mean
end
verifyError(tc, @() swerling_pdf(1, 1, 5), 'swerling_pdf:model');
end

function testManyScatterersGiveExponential(tc)
rng(51); n = 40; pos = 5*(rand(n, 3) - 0.5); pos(:,3) = 0;
phi = (0:0.01:360)';
u = db2lin(rcs_scatterers(zeros(n,1), pos, [1e3*ones(size(phi)), 90*ones(size(phi)), phi], 3e9));
u = u/mean(u);
verifyEqual(tc, std(u), 1, 'AbsTol', 0.1);                               % exponential: std = mean
end

function testSwerlingDetection(tc)
% Pd formulas: limits and Monte Carlo
verifyEqual(tc, pd_swerling3(-80, 1e-4), 1e-4, 'RelTol', 1e-3);          % S -> 0: Pd -> Pfa
verifyGreaterThan(tc, pd_swerling3(20, 1e-6), pd_swerling1(20, 1e-6));   % dominant scatterer helps
rng(52); n = 4e5; S = db2lin(15); Pfa = 1e-4; Y = -log(Pfa);
u = sum(-log(rand(n, 2)), 2)/2;                                          % chi-square 4, unit mean
z = abs(sqrt(S*u).*exp(1i*2*pi*rand(n,1)) + sqrt(1/2)*(randn(n,1) + 1i*randn(n,1))).^2;
verifyEqual(tc, mean(z > Y), pd_swerling3(15, Pfa), 'AbsTol', 0.004);
end
