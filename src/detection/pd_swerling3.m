function Pd = pd_swerling3(SNR_dB, Pfa)
%PD_SWERLING3  Single-pulse Pd of a Swerling III/IV target with a known-noise threshold.  [2026]
%   Pd = pd_swerling3(SNR_dB, Pfa)
%
%   With an RCS distributed as chi-square with 4 dof (mean SNR S) and
%   complex Gaussian noise, the moment-generating function of the
%   square-law output z factorises as (1-t)/(1-beta t)^2, beta = 1 + S/2,
%   i.e. a mixture of an exponential and a Gamma(2) variable. Its tail gives
%       Pd = exp(-Y/beta) * [1 + (1 - 1/beta) Y / beta],   Y = -ln(Pfa)
%   S -> 0 recovers Pd = Pfa. (Swerling I, for comparison: Pd = Pfa^(1/(1+S)).)
S = db2lin(SNR_dB);
Y = -log(Pfa);
beta = 1 + S/2;
Pd = exp(-Y./beta) .* (1 + (1 - 1./beta) .* Y ./ beta);
end
