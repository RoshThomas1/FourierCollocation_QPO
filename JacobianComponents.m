
function [d2u,d_delDc,d_delDs,d_delDcdA,d_delDcdB,d_delDsdA,d_delDsdB] = ...
    JacobianComponents(w,A,B,t,Np,Dc,Ds,dDc,dDs,d2Dc,d2Ds)
% JacobianComponents  Build the 3x3 Jacobian for the NAFF Newton step.
%
%   Computes partial derivatives of the NAFF objective function u(w,A,B)
%   and its cross-terms delta_Dc and delta_Ds w.r.t. [w, A, B].
%   These form the 3x3 system solved at each Newton iteration (Eq. 33).

% Derivatives of delta_D w.r.t. amplitudes A and B
d_delDcdA = 1/Np * cos(w*t) * [cos(w*t)].';
d_delDcdB = 1/Np * sin(w*t) * [cos(w*t)].';
d_delDsdA = 1/Np * cos(w*t) * [sin(w*t)].';
d_delDsdB = 1/Np * sin(w*t) * [sin(w*t)].';

% Derivative of delta_D w.r.t. frequency w
d_delDc = 1/Np * sum((-t*A.*sin(w*t) + t*B.*cos(w*t)).*cos(w*t) + ...
                       (A*cos(w*t)   + B*sin(w*t)).*(-t.*sin(w*t))) - dDc;
d_delDs = 1/Np * sum((-t*A.*sin(w*t) + t*B.*cos(w*t)).*sin(w*t) + ...
                       (A*cos(w*t)   + B*sin(w*t)).*(t.*cos(w*t)))  - dDs;

% Second derivative of NAFF objective u w.r.t. w
d2u = -1/(Dc^2+Ds^2)^(3/2) * (2*Dc*dDc + 2*Ds*dDs) * (Dc*dDc + Ds*dDs) + ...
       2/sqrt(Dc^2+Ds^2)    * (dDc^2 + Dc*d2Dc + dDs^2 + Ds*d2Ds);
end

