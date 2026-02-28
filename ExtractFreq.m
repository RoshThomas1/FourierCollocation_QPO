
function [w,A,B] = ExtractFreq(y_p,t,w_s,Np)
% ExtractFreq  Extract dominant frequency from a windowed signal via FFT.
%
%   Identifies the peak in the FFT magnitude spectrum, then computes the
%   cosine and sine amplitudes at that frequency via inner products.
%
%   INPUTS:
%     y_p - windowed signal samples
%     t   - time vector
%     w_s - frequency axis corresponding to FFT output [rad/TU]
%     Np  - number of samples
%
%   OUTPUT:
%     w   - dominant frequency [rad/TU]
%     A   - cosine amplitude
%     B   - sine amplitude

Y = abs(fftshift(fft(y_p, Np+1)));
w = max(w_s(max(Y) == Y));          % pick positive-frequency peak
A = 1/Np * y_p * [cos(w*t)].';
B = 1/Np * y_p * [sin(w*t)].';
end

