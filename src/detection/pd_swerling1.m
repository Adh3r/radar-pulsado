function Pd = pd_swerling1(SNR_dB, Pfa, M)
%PD_SWERLING1  Single-pulse Pd of a Swerling I/II target, fixed threshold or CA-CFAR.  [2026]
%   Pd = pd_swerling1(SNR_dB, Pfa)      known noise power (fixed threshold):
%                                       Pd = Pfa^(1/(1+SNR))
%   Pd = pd_swerling1(SNR_dB, Pfa, M)   cell-averaging CFAR with M training cells:
%                                       Pd = (1 + alpha/(M (1+SNR)))^(-M)
%   A Swerling I echo plus noise is complex Gaussian with power (1+SNR)
%   times the noise power, so its square-law sample is exponential: the
%   CA-CFAR expression is eq. 6.1 with the noise mean scaled by (1+SNR).
%   M = Inf gives the fixed-threshold result; the SNR gap between the two
%   curves at a given Pd is the CFAR loss.
if nargin < 3, M = Inf; end
snr = db2lin(SNR_dB);
if isinf(M)
    Pd = Pfa.^(1 ./ (1 + snr));
else
    alpha = ca_cfar_alpha(Pfa, M);
    Pd = (1 + alpha ./ (M*(1 + snr))).^(-M);
end
end
