
function y = ReconSig(A,B,w,t)
% ReconSig  Reconstruct a signal from NAFF amplitude/frequency output.

w = w.';
y = A*cos(w*t) + B*sin(w*t);
end

