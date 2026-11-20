function codes = barker_search(K, tol)
%BARKER_SEARCH  Exhaustive search of binary codes of length K with Barker sidelobes (Lab 4.7).
%   codes = barker_search(K)       rows of +1/-1 whose aperiodic autocorrelation
%   codes = barker_search(K, tol)  has every sidelobe |r(k)| <= 1 + tol (default 1e-9)
%   Sequences are enumerated with the first chip fixed to +1 (a code and its
%   negative have the same autocorrelation), i.e. 2^(K-1) candidates.
%   Positive controls: K = 7 and 13 return their Barker codes; K = 6 returns
%   an empty matrix (no Barker code of length 6 exists).
if nargin < 2, tol = 1e-9; end
n = 2^(K-1);
bits  = dec2bin(0:n-1, K-1) - '0';
cands = [ones(n,1), 1 - 2*bits];         % +1 / -1, first chip fixed
keep  = false(n,1);
for i = 1:n
    r = conv(cands(i,:), flip(cands(i,:)));
    r(K) = 0;                             % remove the zero-lag peak
    keep(i) = max(abs(r)) <= 1 + tol;
end
codes = cands(keep, :);
end
