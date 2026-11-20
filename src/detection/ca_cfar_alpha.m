function alpha = ca_cfar_alpha(Pfa, M, N)
%CA_CFAR_ALPHA  CA-CFAR scaling factor for a desired false-alarm probability (Lab 6.2).
%   alpha = ca_cfar_alpha(Pfa, M)      single pulse, closed form:
%                                      alpha = M (Pfa^(-1/M) - 1)
%   alpha = ca_cfar_alpha(Pfa, M, N)   [2026] cells integrating N pulses,
%                                      solved numerically from CA_CFAR_PFA
%   As M -> inf, alpha -> -ln(Pfa): the fixed threshold for known noise power
%   (THRESHOLD_FROM_PFA with unit noise). A finite M needs a larger alpha,
%   because the noise estimate itself is noisy: that excess is the CFAR loss.
if nargin < 3, N = 1; end
alpha = zeros(size(Pfa));
for i = 1:numel(Pfa)
    if N == 1
        alpha(i) = M*(Pfa(i)^(-1/M) - 1);
    else
        f = @(la) log(ca_cfar_pfa(exp(la), M, N)) - log(Pfa(i));   % solve in log(alpha)
        % starting bracket: known-noise threshold relative to the cell mean (N)
        % (the CA-CFAR alpha is always larger; widen the bracket until it is enclosed)
        a0 = threshold_from_pfa(Pfa(i), 1, N) / N;
        lo = log(a0) - 1; hi = log(a0) + 1;
        while f(hi) > 0, hi = hi + 1; end
        alpha(i) = exp(fzero(f, [lo hi]));
    end
end
end
