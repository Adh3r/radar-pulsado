function [y, h] = mti_canceller(x, order, dim)
%MTI_CANCELLER  Single / double (n-th order) delay-line canceller along slow time (Lab 7).
%   [y, h] = mti_canceller(x, order)        order 1: single canceller  y[n] = x[n] - x[n-1]
%   [y, h] = mti_canceller(x, order, dim)   order 2: double canceller  y[n] = x[n] - 2x[n-1] + x[n-2]
%     x     pulse-to-pulse samples: a vector, or a matrix fast time x slow time
%           (default for matrices: dim = 2, each range cell filtered across pulses)
%   Coefficients are binomial with alternating sign, h = (1 - z^-1)^order.
%   Only fully-formed outputs are returned ('valid'): the length along dim
%   drops by 'order'. Output noise power = sum(h.^2) x input noise power
%   (2 for the single, 6 for the double canceller).
k = 0:order;
h = (-1).^k .* arrayfun(@(i) nchoosek(order, i), k);
if nargin < 3
    if isvector(x) && iscolumn(x), dim = 1; else, dim = 2; end
end
if dim == 2
    y = conv2(x, h, 'valid');            % row kernel: filter along dim 2
else
    y = conv2(x, h(:), 'valid');         % column kernel: filter along dim 1
end
end
