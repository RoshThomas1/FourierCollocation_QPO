
function f = QPO_fun_first(x,t_p,M,J,u)
% QPO_fun_first  Objective function for the initial QPO refinement.
%
%   Computes the Jacobi constant residual C(t) - C(t=0) over 1000 points
%   sampled from the Fourier trajectory. A dynamically consistent torus
%   has C = constant everywhere, so this residual should be driven to zero.
%
%   DESIGN VECTOR x = [Ax(M), Ay(M), Az(M), Bx(M), By(M), Bz(M), w0, wP, phi]
%
%   INPUTS:
%     x   - design vector (Fourier coefficients + frequencies + phase)
%     t_p - [1 x N_arc] arc start times
%     M   - number of frequency combinations (length of J)
%     J   - [M x 2] integer combination matrix
%     u   - CR3BP mass parameter
%
%   OUTPUT:
%     f   - [1000 x 1] Jacobi constant residuals C(t) - C(0)

Ax = x(1:M);       Ay = x(M+1:2*M);     Az = x(2*M+1:3*M);
Bx = x(3*M+1:4*M); By = x(4*M+1:5*M);   Bz = x(5*M+1:6*M);
w  = x(end-2:end-1).';

t = linspace(t_p(1), t_p(end)*2, 1000);

% Fourier trajectory evaluation
X  = Ax*cos(J*w*t) + Bx*sin(J*w*t);
Y  = Ay*cos(J*w*t) + By*sin(J*w*t);
Z  = Az*cos(J*w*t) + Bz*sin(J*w*t);
Xd = -(Ax.*(J*w).')*sin(J*w*t) + (Bx.*(J*w).')*cos(J*w*t);
Yd = -(Ay.*(J*w).')*sin(J*w*t) + (By.*(J*w).')*cos(J*w*t);
Zd = -(Az.*(J*w).')*sin(J*w*t) + (Bz.*(J*w).')*cos(J*w*t);

IC_Q = [X(1), Y(1), Z(1), Xd(1), Yd(1), Zd(1)];

% Jacobi constant at each point along the Fourier trajectory
r1 = sqrt((X+u).^2   + Y.^2 + Z.^2);
r2 = sqrt((X-1+u).^2 + Y.^2 + Z.^2);
O  = (1-u)./r1 + u./r2 + 0.5*(X.^2 + Y.^2);
Cq = -(Xd.^2 + Yd.^2 + Zd.^2) + 2*O;

C = JacobiConstant(IC_Q, u);

f = (Cq - C);
end

