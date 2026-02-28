
function y = ReduceSignal(y,t,A,B,w)
% ReduceSignal  Remove one frequency component from a signal.
%
%   Subtracts A*cos(w*t) + B*sin(w*t) from y, leaving the residual for
%   subsequent NAFF iterations.

y = y - (A*cos(w*t) + B*sin(w*t));
end

