function [NF_dB, F, contrib] = cascade_noise_figure(G_dB, NF_dB_stages)
%CASCADE_NOISE_FIGURE  Overall noise figure of cascaded stages (Friis formula).
%   [NF_dB, F, contrib] = cascade_noise_figure(G_dB, NF_dB_stages)
%     G_dB          gain of each stage [dB] (passive stage: G = -loss)
%     NF_dB_stages  noise figure of each stage [dB] (passive stage: NF = loss)
%   Returns the total NF [dB], the total noise factor F and the contribution
%   of each stage referred to the input:
%       F = F1 + (F2-1)/G1 + (F3-1)/(G1 G2) + ...
%   The first stage dominates: its noise is not divided by any gain.
if numel(G_dB) ~= numel(NF_dB_stages)
    error('cascade_noise_figure:size', 'G_dB and NF_dB_stages must have the same length.');
end
G  = db2lin(G_dB(:).');
Fi = db2lin(NF_dB_stages(:).');
Gacc    = [1, cumprod(G(1:end-1))];          % gain in front of each stage
contrib = [Fi(1), (Fi(2:end) - 1) ./ Gacc(2:end)];
F       = sum(contrib);
NF_dB   = lin2db(F);
end
