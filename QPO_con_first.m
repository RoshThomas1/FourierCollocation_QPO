
function [c,ceq] = QPO_con_first(x,t_p,M,J,u,ICp,CM_vec,ds,n_points)
% QPO_con_first  Constraint function for the initial QPO refinement.
%
%   Enforces two sets of constraints:
%
%   INEQUALITY (c):
%     None (c = [])
%
%   EQUALITY (ceq):
%     (1) Center manifold displacement (at t_p(1) only):
%           dot(IC_Q - ICp, xi) = ds
%         where xi = Re(CM_vec)*cos(phi) + Im(CM_vec)*sin(phi)
%         This places the first QPO at distance ds from the PO along the
%         center manifold direction (Eq. 10 of the paper).
%
%     (2) Collocation on each arc [t_p(ii), t_p(ii)+dt]:
%           [X,Y,Z,Xd,Yd,Zd]_Fourier = [x,y,z,xd,yd,zd]_ode113
%         at n_points sample times per arc.
%
%   INPUTS:
%     x        - design vector
%     t_p      - [1 x N_arc] arc start times
%     M        - number of frequency combinations
%     J        - [M x 2] integer combination matrix
%     u        - CR3BP mass parameter
%     ICp      - [1 x 6] reference PO initial condition
%     CM_vec   - [1 x 6] complex center manifold eigenvector
%     ds       - center manifold step size [LU]
%     n_points - collocation points per arc
%
%   OUTPUT:
%     c    - inequality constraints (empty)
%     ceq  - equality constraints

Ax = x(1:M);       Ay = x(M+1:2*M);     Az = x(2*M+1:3*M);
Bx = x(3*M+1:4*M); By = x(4*M+1:5*M);   Bz = x(5*M+1:6*M);
w  = x(end-2:end-1).'; phi = x(end);
dt = t_p(2) - t_p(1);
ceq = [];

for ii = 1:length(t_p)
    t = linspace(t_p(ii), t_p(ii)+dt, n_points);

    X  = Ax*cos(J*w*t) + Bx*sin(J*w*t);
    Y  = Ay*cos(J*w*t) + By*sin(J*w*t);
    Z  = Az*cos(J*w*t) + Bz*sin(J*w*t);
    Xd = -(Ax.*(J*w).')*sin(J*w*t) + (Bx.*(J*w).')*cos(J*w*t);
    Yd = -(Ay.*(J*w).')*sin(J*w*t) + (By.*(J*w).')*cos(J*w*t);
    Zd = -(Az.*(J*w).')*sin(J*w*t) + (Bz.*(J*w).')*cos(J*w*t);

    % Center manifold direction rotated by phase angle phi
    xi = real(CM_vec)*cos(phi) + imag(CM_vec)*sin(phi);
    xi = xi / norm(xi);

    IC_Q = [X(1), Y(1), Z(1), Xd(1), Yd(1), Zd(1)];

    if ii == 1
        % Enforce prescribed center manifold displacement
        r    = dot([IC_Q - ICp], xi);
        c    = [];
        ceq  = [100*(r - ds)];
    end

    options = odeset('RelTol',1e-8,'AbsTol',1e-8);
    [t_q,s] = ode113(@(t,s) CR3BP(t,s,u), t, IC_Q, options);

    if t_q(end) == t(end)
        % Collocation: Fourier trajectory must match CR3BP propagation
        ceq = [ceq, (X-s(:,1).'), (Y-s(:,2).'), (Z-s(:,3).'), ...
                     Xd-s(:,4).', Yd-s(:,5).', Zd-s(:,6).'];
    else
        continue;   % skip arc if integrator failed to reach end time
    end
end
end

