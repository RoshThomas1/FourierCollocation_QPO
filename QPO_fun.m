
function f = QPO_fun(x,t_p,M,J,u)
% QPO_fun  Objective function for QPO continuation and harmonic addition.
%
%   Identical to QPO_fun_first but without the phase variable phi in the
%   design vector (removed after the first QPO is found).
%
%   NOTE: Xdd, Ydd, Zdd are computed below but not used in the residual.
%   These are second derivatives of the Fourier trajectory, retained from
%   an earlier formulation. They do not affect the result and may be
%   removed in a future cleanup.
%
%   INPUTS / OUTPUT: see QPO_fun_first (phi removed from x)

Ax = x(1:M);       Ay = x(M+1:2*M);     Az = x(2*M+1:3*M);
Bx = x(3*M+1:4*M); By = x(4*M+1:5*M);   Bz = x(5*M+1:6*M);
w  = x(end-1:end).';

t = linspace(t_p(1), t_p(end)*2, 1000);

X  = Ax*cos(J*w*t) + Bx*sin(J*w*t);
Y  = Ay*cos(J*w*t) + By*sin(J*w*t);
Z  = Az*cos(J*w*t) + Bz*sin(J*w*t);
Xd = -(Ax.*(J*w).')*sin(J*w*t) + (Bx.*(J*w).')*cos(J*w*t);
Yd = -(Ay.*(J*w).')*sin(J*w*t) + (By.*(J*w).')*cos(J*w*t);
Zd = -(Az.*(J*w).')*sin(J*w*t) + (Bz.*(J*w).')*cos(J*w*t);

IC_Q = [X(1), Y(1), Z(1), Xd(1), Yd(1), Zd(1)];

% NOTE: Second derivatives computed but unused in residual (dead code).
Xdd = -(Ax.*((J*w).^2).')*cos(J*w*t) - (Bx.*((J*w).^2).')*sin(J*w*t);
Ydd = -(Ay.*((J*w).^2).')*cos(J*w*t) - (By.*((J*w).^2).')*sin(J*w*t);
Zdd = -(Az.*((J*w).^2).')*cos(J*w*t) - (Bz.*((J*w).^2).')*sin(J*w*t);

r1 = sqrt((X+u).^2   + Y.^2 + Z.^2);
r2 = sqrt((X-1+u).^2 + Y.^2 + Z.^2);
O  = (1-u)./r1 + u./r2 + 0.5*(X.^2 + Y.^2);
Cq = -(Xd.^2 + Yd.^2 + Zd.^2) + 2*O;

C = JacobiConstant(IC_Q, u);
f = (Cq - C);
end

