function Q = marcum_q(a, b)
%MARCUM_Q  First-order Marcum Q function Q1(a, b), base MATLAB only.  [2026]
%   Q = marcum_q(a, b)   (a, b >= 0, arrays broadcast)
%
%       Q1(a,b) = integral_b^inf  x exp(-(x^2 + a^2)/2) I0(a x) dx
%
%   Evaluated directly on the upper tail with the exponentially scaled
%   Bessel function, so tiny values (Pfa ~ 1e-12) keep full relative accuracy.
%   Q1(0, b) = exp(-b^2/2) is the Rayleigh tail.
[a, b] = deal(a + 0*b, b + 0*a);                 % broadcast to a common size
Q = zeros(size(a));
for i = 1:numel(a)
    ai = a(i);
    f = @(x) x .* exp(-(x - ai).^2/2) .* besseli(0, ai*x, 1);   % I0(z) = besseli(0,z,1)*e^z
    Q(i) = integral(f, b(i), Inf, 'RelTol', 1e-10, 'AbsTol', 1e-300);
end
Q = min(Q, 1);
end
