
function [Dc,Ds,dDc,dDs,d2Dc,d2Ds] = ShootingParams(y,t,w,Np)
% ShootingParams  Compute inner products needed for NAFF Newton iteration.
%
%   Evaluates Dc, Ds and their first and second derivatives with respect
%   to frequency w, used to solve the NAFF optimization problem (Eq. 32).
%
%   All quantities are normalized by 1/Np (discrete inner product).

dDc  = 1/Np * y * [-sin(w*t).*t].';
dDs  = 1/Np * y * [ cos(w*t).*t].';
Dc   = 1/Np * y * [ cos(w*t)   ].';
Ds   = 1/Np * y * [ sin(w*t)   ].';
d2Dc = 1/Np * y * [-cos(w*t).*t.^2].';
d2Ds = 1/Np * y * [-sin(w*t).*t.^2].';
end

