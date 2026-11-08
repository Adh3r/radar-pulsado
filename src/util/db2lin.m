function y = db2lin(x_dB)
%DB2LIN  Power ratio in dB -> linear (base-MATLAB replacement for db2pow).
y = 10.^(x_dB/10);
end
