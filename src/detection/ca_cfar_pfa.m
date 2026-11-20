function Pfa = ca_cfar_pfa(alpha, M, N)
%CA_CFAR_PFA  False-alarm probability of a cell-averaging CFAR (Lab 6, eq. 6.1).
%   Pfa = ca_cfar_pfa(alpha, M)      single pulse:  Pfa = (1 + alpha/M)^(-M)
%   Pfa = ca_cfar_pfa(alpha, M, N)   [2026] each cell is the sum of N pulses
%                                    (non-coherent integration)
%     alpha  threshold scaling factor (T = alpha * mean of the training cells)
%     M      number of training cells
%
%   Square-law samples of complex Gaussian noise are exponential, so the CUT
%   is Gamma(N) and the sum of the M training cells is Gamma(M N). Averaging
%   over the training sum gives, with beta = alpha/M,
%       Pfa = sum_{k=0}^{N-1} C(MN+k-1, k) beta^k / (1 + beta)^(MN+k)
%   which does not depend on the noise power: that is the CFAR property.
%   N = 1 reduces to eq. 6.1 of the lab.
if nargin < 3, N = 1; end
beta = alpha ./ M;
Pfa = zeros(size(beta));
for k = 0:N-1
    logterm = gammaln(M*N + k) - gammaln(k + 1) - gammaln(M*N) ...
              + k*log(beta) - (M*N + k)*log1p(beta);
    Pfa = Pfa + exp(logterm);
end
end
