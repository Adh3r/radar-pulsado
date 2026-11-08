function y_dB = lin2db(x)
%LIN2DB  Linear power ratio -> dB (base-MATLAB replacement for pow2db).
y_dB = 10*log10(x);
end
