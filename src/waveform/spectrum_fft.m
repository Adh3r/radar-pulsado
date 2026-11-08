function [f, X] = spectrum_fft(x, fs)
%SPECTRUM_FFT  Two-sided magnitude spectrum normalised by N (fft + fftshift).
%   [f, X] = spectrum_fft(x, fs) -> frequency axis f [Hz] centred on 0, and |X|.
%   Dividing by N makes the height of spectral lines independent of the
%   record length; fftshift moves the negative frequencies to the left.
N = numel(x);
X = abs(fftshift(fft(x)))/N;
f = (-floor(N/2):ceil(N/2)-1) * fs/N;
end
