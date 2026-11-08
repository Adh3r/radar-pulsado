function [SNR_dB, SNR_lin] = radar_snr(Pt, G_dBi, fc, B, NF_dB, L_dB, sigma, R)
%RADAR_SNR  Single-pulse SNR of a monostatic radar (radar range equation).
%   [SNR_dB, SNR_lin] = radar_snr(Pt, G_dBi, fc, B, NF_dB, L_dB, sigma, R)
%     Pt peak power [W], G_dBi antenna gain (Tx = Rx), fc carrier [Hz],
%     B noise bandwidth [Hz], NF_dB receiver noise figure, L_dB extra losses,
%     sigma RCS [m^2], R target range(s) [m]. Inputs broadcast (implicit expansion).
%
%              Pt G^2 lambda^2 sigma
%      SNR = -------------------------
%            (4 pi)^3 R^4 k T0 B F L
%
%   Two spherical spreadings (out and back) give R^4: doubling range costs 12 dB.
k = radar_constants();
lambda = k.c./fc;
Pr = Pt.*db2lin(G_dBi).^2.*lambda.^2.*sigma ./ ((4*pi)^3 * R.^4 .* db2lin(L_dB));
N  = k.kB*k.T0.*B.*db2lin(NF_dB);
SNR_lin = Pr./N;
SNR_dB  = lin2db(SNR_lin);
end
