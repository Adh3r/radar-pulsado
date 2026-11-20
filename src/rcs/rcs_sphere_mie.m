function s_norm = rcs_sphere_mie(ka)
%RCS_SPHERE_MIE  Normalised monostatic RCS of a perfectly conducting sphere (Mie series).
%   s_norm = rcs_sphere_mie(ka)      sigma / (pi a^2) for electrical size ka = 2 pi a / lambda
%
%   Exact series for a PEC sphere (e.g. Knott, Ruck):
%     sigma / (pi a^2) = | (1/ka) sum_{n>=1} (-1)^n (2n+1) (b_n - a_n) |^2
%     a_n = j_n(ka) / h_n(ka)
%     b_n = [ka j_{n-1}(ka) - n j_n(ka)] / [ka h_{n-1}(ka) - n h_n(ka)]
%   with spherical Bessel j_n and Hankel h_n = j_n + i y_n, built from MATLAB's
%   besselj / bessely (J_{n+1/2}). Truncated at N = ka + 4 ka^(1/3) + 3 terms.
%   Limits: Rayleigh region sigma/(pi a^2) -> 9 (ka)^4, optical region -> 1.
s_norm = zeros(size(ka));
for i = 1:numel(ka)
    x = ka(i);
    N = ceil(x + 4*x^(1/3) + 3);
    n = 0:N;
    jn = sqrt(pi/(2*x)) * besselj(n + 0.5, x);
    yn = sqrt(pi/(2*x)) * bessely(n + 0.5, x);
    hn = jn + 1i*yn;
    m  = 1:N;                                    % orders n >= 1 (index m+1 holds order m)
    a = jn(m+1) ./ hn(m+1);
    b = (x*jn(m) - m.*jn(m+1)) ./ (x*hn(m) - m.*hn(m+1));
    s_norm(i) = abs(sum((-1).^m .* (2*m + 1) .* (b - a)) / x)^2;
end
end
