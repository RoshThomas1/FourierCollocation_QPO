

function C = JacobiConstant(s,u)
% JacobiConstant  Compute the Jacobi constant for CR3BP states.
%
%   C = 2*Omega - v^2, where Omega is the effective potential.
%   The Jacobi constant is the sole integral of motion in the CR3BP.
%
%   INPUTS:
%     s - [N x 6] matrix of states [x, y, z, xdot, ydot, zdot]
%     u - CR3BP mass parameter mu
%
%   OUTPUT:
%     C - [N x 1] Jacobi constant values

r1 = sqrt((s(:,1)+u).^2    + s(:,2).^2 + s(:,3).^2);
r2 = sqrt((s(:,1)-1+u).^2  + s(:,2).^2 + s(:,3).^2);
O  = (1-u)./r1 + u./r2 + 0.5*(s(:,1).^2 + s(:,2).^2);
C  = -(( s(:,4).^2 + s(:,5).^2 + s(:,6).^2 ) - 2*O);
end
