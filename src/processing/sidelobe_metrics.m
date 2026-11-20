function m = sidelobe_metrics(y)
%SIDELOBE_METRICS  Peak/sidelobe figures of a compressed pulse (Lab 4.3 metrics).
%   m = sidelobe_metrics(y)    y = matched-filter output (complex or real)
%   Fields:
%     peak        main-lobe peak |y|
%     PSL_dB      peak sidelobe level, 20 log10(max sidelobe / peak)
%     ISL_dB      integrated sidelobe level, 10 log10(sum sidelobe^2 / peak^2)
%     n_lobes     number of lobes (local maxima above 1e-6 * peak, main included)
%     mainlobe_width  null-to-null width of the main lobe [samples]
%     length      output length [samples]
a = abs(y(:));
[m.peak, ip] = max(a);
m.length = numel(a);

% main lobe: walk downhill from the peak on both sides
il = ip; while il > 1        && a(il-1) < a(il), il = il - 1; end
ir = ip; while ir < numel(a) && a(ir+1) < a(ir), ir = ir + 1; end
m.mainlobe_width = ir - il;

side = a; side(il:ir) = 0;
m.PSL_dB = 20*log10(max(side) / m.peak);
m.ISL_dB = 10*log10(sum(side.^2) / m.peak^2);

isMax = [a(1) > a(2); a(2:end-1) > a(1:end-2) & a(2:end-1) >= a(3:end); a(end) > a(end-1)];
m.n_lobes = nnz(isMax & a > 1e-6*m.peak);
end
