function Pd = pd_swerling0(SNR_dB, Pfa)
%PD_SWERLING0  Exact single-pulse Pd of a non-fluctuating target (envelope detector).  [2026]
%   Pd = pd_swerling0(SNR_dB, Pfa)   (arrays broadcast)
%   Complex Gaussian noise of power sigma^2 and a constant-amplitude echo A
%   (SNR = A^2/sigma^2). Normalising by the per-quadrature std sigma/sqrt(2):
%       Pd = Q1( sqrt(2 SNR), sqrt(-2 ln Pfa) )
%   This is the curve that Albersheim's equation approximates.
Pd = marcum_q(sqrt(2*db2lin(SNR_dB)), sqrt(-2*log(Pfa)));
end
