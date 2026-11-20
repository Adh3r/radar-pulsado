function [Pfa, Pd] = roc_empirical(h1, h0, T)
%ROC_EMPIRICAL  Empirical ROC from detector outputs under H1 and H0.  [2026]
%   [Pfa, Pd] = roc_empirical(h1, h0, T)
%     h1  detector statistics with target present (any array)
%     h0  detector statistics with noise only (any array)
%     T   thresholds (vector)
%   Pd(k) = fraction of h1 above T(k); Pfa(k) = fraction of h0 above T(k).
Pd  = zeros(size(T));
Pfa = zeros(size(T));
for k = 1:numel(T)
    Pd(k)  = mean(h1(:) > T(k));
    Pfa(k) = mean(h0(:) > T(k));
end
end
