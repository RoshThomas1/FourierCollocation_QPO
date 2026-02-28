
function [c,ceq] = QPO_con_added(x,t_p,M,J,u,IC_pert,n_points)
% QPO_con_added  Constraint function for harmonic addition refinement.
%
%   Used when re-solving the collocation problem after new harmonics have
%   been appended to the frequency basis. Unlike QPO_con, no arclength
%   or phase constraints are applied -- instead, the first arc's initial
%   condition is pinned to a known good initial condition IC_pert from
%   the original continuation family, anchoring the solution in phase space.
%
%   EQUALITY (ceq):
%     (1) Initial condition pinning (at ii=1):
%           IC_Q = IC_pert   (scaled by 50 for numerical conditioning)
%     (2) Collocation on each arc (same as other constraint functions)
%
%   INPUTS:
%     x        - design vector
%     t_p      - arc start times
%     M        - number of frequency combinations
%     J        - [M x 2] integer combination matrix
%     u        - CR3BP mass parameter
%     IC_pert  - [1 x 6] known initial condition to anchor the solution
%     n_points - collocation points per arc

Ax = x(1:M);       Ay = x(M+1:2*M);     Az = x(2*M+1:3*M);
Bx = x(3*M+1:4*M); By = x(4*M+1:5*M);   Bz = x(5*M+1:6*M);
w  = x(end-1:end).';
E  = Ax*Ax.' + Ay*Ay.' + Az*Az.' + Bx*Bx.' + By*By.' + Bz*Bz.';
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

    IC_Q = [X(1), Y(1), Z(1), Xd(1), Yd(1), Zd(1)];

    if ii == 1
        c   = [];
        % Pin initial condition to known value from continuation family
        ceq = 50 * [IC_pert - IC_Q];
    end

    options = odeset('RelTol',1e-8,'AbsTol',1e-8);
    [t_q,s] = ode113(@(t,s) CR3BP(t,s,u), t, IC_Q, options);

    if t_q(end) == t(end)
        ceq = [ceq, (X-s(:,1).'), (Y-s(:,2).'), (Z-s(:,3).'), ...
                     Xd-s(:,4).', Yd-s(:,5).', Zd-s(:,6).'];
    else
        continue;
    end
end
end
