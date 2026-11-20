function [idx, val] = detect_peaks(p, T)
%DETECT_PEAKS  Threshold detection with one report per contiguous exceedance.  [2026]
%   [idx, val] = detect_peaks(p, T)
%     p  detector output along range (vector, e.g. |y|^2)
%     T  threshold (scalar or same size as p)
%   Each run of consecutive samples above T is reported once, at its maximum,
%   so a compressed pulse spread over a few samples gives a single detection.
p = p(:); above = p > T(:);
d = diff([0; above; 0]);
starts = find(d == 1); stops = find(d == -1) - 1;
idx = zeros(numel(starts), 1); val = idx;
for k = 1:numel(starts)
    [val(k), j] = max(p(starts(k):stops(k)));
    idx(k) = starts(k) + j - 1;
end
end
