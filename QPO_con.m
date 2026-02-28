
function [c,ceq] = QPO_con(x,t_p,M,J,u,ds,coeff_xi_old,x_old,n_points)
% QPO_con  Constraint function for QPO arclength continuation.
%
%   Enforces three constraint types (evaluated only at ii=1 for the
%   algebraic constraints, and on each arc for the collocation):
%
%   INEQUALITY (c):
%     Fourier energy monotonicity: norm(F_new)^2 >= norm(F_old)^2
%     Penalized as: 1e5 * max(0, E_old - E_new)
%     This enforces that QPO amplitude increases along the family.
%
%   EQUALITY (ceq):
%     (1) Arclength constraint (Eq. 13):
%           (x - x_old) . coeff_xi_old = ds
%         Steps a fixed arclength ds in Fourier coefficient space.
%
%     (2,3) Phase conditions (Eq. 14-15):
%           (F - F_old) . dQ1 = 0
%           (F - F_old) . dQ2 = 0
%         Prevents the continuation from drifting to a phase-equivalent
%         torus. dQ1 and dQ2 are Fourier-space representations of the
%         partial derivatives d/dQ1 and d/dQ2 of the previous solution.
%
%     (4) Collocation on each arc (same as QPO_con_first).
%
%   INPUTS:
%     x            - current design vector
%     t_p          - arc start times
%     M            - number of frequency combinations
%     J            - [M x 2] integer combination matrix
%     u            - CR3BP mass parameter
%     ds           - arclength step size
%     coeff_xi_old - previous arclength direction (weighted)
%     x_old        - previous solution vector
%     n_points     - collocation points per arc

Ax = x(1:M);       Ay = x(M+1:2*M);     Az = x(2*M+1:3*M);
Bx = x(3*M+1:4*M); By = x(4*M+1:5*M);   Bz = x(5*M+1:6*M);
w  = x(end-1:end).';
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

    % Retrieve previous Fourier coefficients for phase condition computation
    Axo = x_old(1:M);       Ayo = x_old(M+1:2*M);     Azo = x_old(2*M+1:3*M);
    Bxo = x_old(3*M+1:4*M); Byo = x_old(4*M+1:5*M);   Bzo = x_old(5*M+1:6*M);

    % Phase condition direction vectors in Fourier coefficient space
    % dQ1 = d(q)/d(Q1): derivative w.r.t. first angle, encoded in coefficients
    dAxo = Axo.*J(:,1).'; dAyo = Ayo.*J(:,1).'; dAzo = Azo.*J(:,1).';
    dBxo = Bxo.*J(:,1).'; dByo = Byo.*J(:,1).'; dBzo = Bzo.*J(:,1).';
    dQ1  = [dBxo, dByo, dBzo, -dAxo, -dAyo, -dAzo]; dQ1 = dQ1/norm(dQ1);

    % dQ2 = d(q)/d(Q2): derivative w.r.t. second angle
    dAxo = Axo.*J(:,2).'; dAyo = Ayo.*J(:,2).'; dAzo = Azo.*J(:,2).';
    dBxo = Bxo.*J(:,2).'; dByo = Byo.*J(:,2).'; dBzo = Bzo.*J(:,2).';
    dQ2  = [dBxo, dByo, dBzo, -dAxo, -dAyo, -dAzo]; dQ2 = dQ2/norm(dQ2);

    F_old = x_old(1:end-2);    % previous Fourier coefficients (excluding w)
    F     = x(1:end-2);        % current  Fourier coefficients (excluding w)

    if ii == 1
        % Arclength constraint
        r = (x - x_old) * coeff_xi_old.';

        % Index excluding DC terms and frequencies (Fourier energy subset)
        ind = 1:1:length(x);
        indices = [[0:2]*length(J)+1, length(x)-[1,0]];
        ind(indices) = [];

        E_old = norm(x_old(ind), 2)^2;
        E     = norm(x(ind),     2)^2;

        % Inequality: Fourier energy must not decrease
        c   = [1e5 * max(0, E_old - E)];

        % Equality: arclength + two phase conditions
        ceq = 100 * [r-ds, dot(F-F_old, dQ1), dot(F-F_old, dQ2)];
    end

    options = odeset('RelTol',1e-8,'AbsTol',1e-8);
    [t_q,s] = ode113(@(t,s) CR3BP(t,s,u), t, IC_Q, options);

    if t_q(end) == t(end)
        ceq = [ceq, X-s(:,1).', Y-s(:,2).', Z-s(:,3).', ...
                    Xd-s(:,4).', Yd-s(:,5).', Zd-s(:,6).'];
    else
        continue;
    end
end
end

