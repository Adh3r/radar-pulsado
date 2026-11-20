function [det, T, noise_est] = ca_cfar(p, n_train, n_guard, alpha)
%CA_CFAR  Cell-averaging CFAR detector on square-law samples (Lab 6.3).
%   [det, T, noise_est] = ca_cfar(p, n_train, n_guard, alpha)
%     p        square-law detector output |y|^2: vector, or matrix of
%              range cells x columns (pulses / trials), processed per column
%     n_train  total number of training cells M (even: M/2 on each side)
%     n_guard  guard cells on EACH side of the cell under test (CUT)
%     alpha    scaling factor, e.g. from CA_CFAR_ALPHA
%   Outputs, same size as p:
%     noise_est  mean of the M training cells around each CUT
%     T          adaptive threshold alpha * noise_est
%     det        logical detections p > T
%   Cells whose window does not fit inside the data (the first and last
%   M/2 + n_guard cells) get NaN and are never declared.
%
%   Window around the CUT (o): [ M/2 training | G guard | o | G guard | M/2 training ]
%   The sliding average is a convolution with that window, i.e. the
%   "sliding-window coefficient matrix" suggested in the lab.
if mod(n_train, 2), error('ca_cfar:train', 'n_train must be even.'); end
h = n_train/2;
w = [ones(h,1); zeros(n_guard,1); 0; zeros(n_guard,1); ones(h,1)] / n_train;
isrowIn = isrow(p);
if isrowIn, p = p(:); end

noise_est = conv2(p, w, 'same');
half = h + n_guard;
noise_est([1:half, end-half+1:end], :) = NaN;     % incomplete windows
T   = alpha * noise_est;
det = p > T;                                       % NaN comparisons are false

if isrowIn, det = det.'; T = T.'; noise_est = noise_est.'; end
end
