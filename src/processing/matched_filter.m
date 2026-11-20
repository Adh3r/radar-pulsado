function y = matched_filter(rx, ref, mode)
%MATCHED_FILTER  Matched-filter (pulse compression) output, h(t) = ref*(-t).
%   y = matched_filter(rx, ref)          full convolution (Lab 4.3 convention)
%   y = matched_filter(rx, ref, 'aligned')
%     rx   received samples: a vector, or a matrix fast-time x pulses
%          (each column filtered independently)
%     ref  transmitted baseband pulse (vector)
%   'full'    : length numel(rx) + numel(ref) - 1; the peak for an echo that
%               starts at sample n0 is at n0 + numel(ref) - 1.
%   'aligned' : same length as rx, delay removed so the peak sits at n0
%               (fast-time sample <-> range bin).
%   The output peak equals the pulse energy sum|ref|^2 and the SNR gain is
%   the time-bandwidth product (numel(ref) for white noise).
if nargin < 3, mode = 'full'; end
h = conj(flip(ref(:)));
L = numel(h);
if isvector(rx)
    y = conv(rx(:), h);
    if isrow(rx), y = y.'; end
else
    y = conv2(rx, h);
end
if strcmpi(mode, 'aligned')
    if isvector(rx)
        y = y(L:L+numel(rx)-1);
    else
        y = y(L:L+size(rx,1)-1, :);
    end
end
end
