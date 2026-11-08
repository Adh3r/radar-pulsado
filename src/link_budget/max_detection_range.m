function R = max_detection_range(Pt, G_dBi, fc, B, NF_dB, L_dB, sigma, SNR_req_dB)
%MAX_DETECTION_RANGE  Range at which the radar equation yields the required SNR.  [2026]
%   R = max_detection_range(Pt, G_dBi, fc, B, NF_dB, L_dB, sigma, SNR_req_dB)
%   Closed-form inversion of RADAR_SNR (SNR ~ R^-4):
%       R = R_ref * 10^((SNR(R_ref) - SNR_req) / 40)
%   Inputs broadcast, so e.g. SNR_req_dB from ALBERSHEIM_SNR can be a matrix.
R_ref   = 1e3;
SNR_ref = radar_snr(Pt, G_dBi, fc, B, NF_dB, L_dB, sigma, R_ref);
R = R_ref * 10.^((SNR_ref - SNR_req_dB)/40);
end
