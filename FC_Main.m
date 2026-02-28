% =========================================================================
% Fourier-Collocation Framework for QPO Family Construction in the CR3BP
% =========================================================================
%
% AUTHORS:  Evan Lantier
%           Department of Aerospace Engineering, Penn State University
%
% DESCRIPTION:
%   Constructs quasi-periodic orbit (QPO) families in the Earth-Moon
%   Circular Restricted Three-Body Problem (CR3BP) using a Fourier-
%   collocation approach. Starting from a reference periodic orbit (PO),
%   the algorithm:
%     1. Differentially corrects the PO to machine precision
%     2. Extracts fundamental frequencies via Floquet theory and NAFF
%     3. Constructs an initial Fourier-based torus approximation
%     4. Refines the torus to enforce CR3BP dynamical consistency
%     5. Continues the QPO family via arclength continuation
%     6. Improves accuracy by systematically adding missing harmonics
%
% UNITS:
%   Length : LU = 384,400 km  (Earth-Moon distance)
%   Time   : TU = 375,197 s   (1/2pi of the synodic period)
%   All states are in normalized (nondimensional) CR3BP units unless
%   explicitly scaled for plotting.
%
% DEPENDENCIES:
%   CR3BP.m, CR3BP_STM.m      -- equations of motion (defined below)
%   planet3D()                -- optional; comment out if not installed
%   CVX + Gurobi              -- required for harmonic addition (Sec 3.5)
%   violinplot()              -- required for error visualization
%
% REFERENCES:
%   Lantier & Eapen, "A Fourier-Collocation Framework for Constructing
%   Quasi-Periodic Orbit Families in the CR3BP," JASS (2024)
%
% =========================================================================

clear all; clc; close all

% =========================================================================
%% SECTION 1: Periodic Orbit Initialization
% =========================================================================
% Define the reference periodic orbit (PO) that serves as the seed for the
% QPO family. The initial condition ICp and period Tp should correspond to
% a known PO in the CR3BP (e.g. from a continuation database or prior run).
%
% The differential corrector below enforces periodicity by targeting
% y(Tp) = 0 and ydot(Tp) = 0, correcting x0, ydot0, and Tp.
% This exploits the z-symmetry of the CR3BP (xz-plane symmetry).
% =========================================================================

% --- Numerical integration tolerances ---
options = odeset('RelTol',1e-10,'AbsTol',1e-10);

% --- Physical scaling constants ---
LU = 384400;                        % km per length unit (Earth-Moon distance)
TU = 3.751971301332566e+05;         % seconds per time unit

% --- Discretization for PO propagation ---
Np = 2^14;                          % number of time nodes along the PO

% --- Initial conditions for the reference PO [x, y, z, xdot, ydot, zdot] ---
% NOTE: This orbit is a vertical orbit near the Earth-Moon L2 point.
%       Modify ICp and Tp to target a different reference orbit.
ICp = [1.153555180492936, 0, 0.139438085327783, 0, -0.215370605225172, 0];
Tp  = 3.221357403048232;            % initial guess for orbital period [TU]

% --- CR3BP mass parameter ---
u = 0.012153593988908;              % mu = m_Moon / (m_Earth + m_Moon)

% -------------------------------------------------------------------------
% Differential Corrector (Newton iteration)
% -------------------------------------------------------------------------
% Corrects [x0, ydot0, Tp] to enforce y(Tp) = 0 and ydot(Tp) = 0.
% The design vector is X = [x0, ydot0, T].
% Convergence is declared when the constraint norm falls below tol.
%
% Uses pseudoinverse of the 2x3 Jacobian DF (underdetermined system),
% which finds the minimum-norm correction at each step.
% -------------------------------------------------------------------------
error = 1; tol = 1e-12;
while error > tol
    t = linspace(0, Tp, Np);
    % Propagate state + State Transition Matrix (STM) simultaneously
    [~,s] = ode45(@(t,s) CR3BP_STM(t,s,u), t, [ICp, reshape(eye(6),1,36)], options);

    % Constraint vector: enforce y and ydot periodicity
    F = [s(end,2), s(end,4)].';

    % Extract terminal state derivative and STM for Jacobian construction
    sdot = CR3BP(0, s(end,1:6), u);
    phi  = reshape(s(end,7:end), 6, 6);    % monodromy matrix at time Tp

    % Build Jacobian DF = d[y(Tp), ydot(Tp)] / d[x0, ydot0, Tp]
    DF = [phi(2,1), phi(2,5), sdot(2);
          phi(4,1), phi(4,5), sdot(4)];

    % Minimum-norm Newton step
    df = pinv(DF) * F;

    % Update design variables
    ICp(1) = ICp(1) - df(1);
    ICp(5) = ICp(5) - df(2);
    Tp     = Tp     - df(3);

    error = norm(F);
    if norm(df) < 1e-11; break; end     % step size convergence fallback
end

% =========================================================================
%% SECTION 2: Monodromy Matrix and Poincare Exponents
% =========================================================================
% Propagate the corrected PO once more to extract the monodromy matrix
% Phi(Tp,0). Its eigenvalues (Floquet multipliers) encode the stability
% of the orbit. The Poincare exponents w_P are defined as:
%
%   w_P = (1/Tp) * log(eigenvalues of Phi)
%
% For a center manifold eigenvalue lambda = exp(i*omega_P*Tp), the
% imaginary part of w_P gives the second fundamental frequency omega_P
% of the surrounding QPO family.
%
% NOTE: The center manifold eigenvectors (columns 2 and 3 below) are
% hardcoded by index. This is valid for the vertical orbit used here,
% where the center pair is consistently the 2nd/3rd eigenvalue. For
% other orbit families, verify the index against the imaginary part of
% w_P before proceeding.
% =========================================================================

[~,s]  = ode45(@(t,s) CR3BP_STM(t,s,u), [0,Tp], [ICp, reshape(eye(6),1,36)], options);
phi    = reshape(s(end,7:end), 6, 6);
[vec,val] = eig(phi);
w_P = 1/Tp * log(diag(val))           % display Poincare exponents

% Extract center manifold eigenvectors and quasi-periodic frequency
% WARNING: indices 2:3 are hardcoded -- verify for non-vertical orbits
vec_center = vec(:,2:3);
w_P        = abs(w_P(3));              % second fundamental frequency [rad/TU]

% Jacobi constant of the reference PO (used later for error tracking)
C_PO = JacobiConstant(ICp, u);

% =========================================================================
%% SECTION 3: Propagate PO and Plot
% =========================================================================
% Propagate the corrected PO over one period at high resolution for
% spectral decomposition and plotting.
% =========================================================================

T = Tp;
t = linspace(0, T, Np);
[~,s] = ode45(@(t,s) CR3BP_STM(t,s,u), t, [ICp, reshape(eye(6),1,36)], options);

r = s(:,1:3).';     % position history [3 x Np]
v = s(:,4:6).';     % velocity history [3 x Np]

% --- Plot the reference periodic orbit ---
% NOTE: planet3D() is a third-party function. Comment out the planet3D
%       calls if it is not installed; the orbit will still plot correctly.
scale = LU / 10^5;  % convert nondimensional LU to units of 10^5 km

figure()
hold on
plot3(scale*r(1,:), scale*r(2,:), scale*r(3,:), 'b', 'LineWidth', 2, 'HandleVisibility', 'off')
opts.Position = scale * [1-u; 0; 0];
opts.Units    = 'EM_km_scale';
planet3D('Moon', opts);
xlabel('x [\times10^5 km]'); ylabel('y [\times10^5 km]'); zlabel('z [\times10^5 km]');
axis equal; view(-45, 15); fontsize(16,'points'); grid on; axis tight;

% =========================================================================
%% SECTION 4: Spectral Decomposition of the Periodic Orbit
% =========================================================================
% Fit a truncated Fourier series to the PO trajectory using only integer
% multiples of the fundamental frequency w0 = 2*pi/Tp. This validates
% that the PO is well-represented by M harmonics before constructing the
% QPO torus.
%
% The least-squares fit uses the normal equations:
%   bin = pinv(X*X') * X * r'
% where X is the [cos; sin] basis matrix evaluated at the sample times.
%
% Residuals are plotted to confirm spectral accuracy.
% =========================================================================

M      = 10;           % maximum harmonic order for PO decomposition
w_fund = 2*pi / Tp;   % fundamental frequency of the PO [rad/TU]

AA = [0:M].';
Jx = AA; Jy = AA; Jz = AA;

% Frequency arrays for each coordinate (all equal for a PO)
w_lin = [{(Jx*w_fund).'}; {(Jy*w_fund).'}; {(Jz*w_fund).'}];

% Fit Fourier coefficients for x
X   = [cos(Jx*w_fund*t); sin(Jx*w_fund*t)];
bin = pinv(X*X.') * X * r(1,:).';
Ax  = bin(1:length(Jx)).';   Bx = bin(length(Jx)+1:end).';

% Fit Fourier coefficients for y
Y   = [cos(Jy*w_fund*t); sin(Jy*w_fund*t)];
bin = pinv(Y*Y.') * Y * r(2,:).';
Ay  = bin(1:length(Jy)).';   By = bin(length(Jy)+1:end).';

% Fit Fourier coefficients for z
Z   = [cos(Jz*w_fund*t); sin(Jz*w_fund*t)];
bin = pinv(Z*Z.') * Z * r(3,:).';
Az  = bin(1:length(Jz)).';   Bz = bin(length(Jz)+1:end).';

A      = [{Ax}; {Ay}; {Az}];
B      = [{Bx}; {By}; {Bz}];
r_lin  = ConstructOrbit(A, B, w_lin, t);

% Plot coordinate residuals between true PO and Fourier reconstruction
figure()
hold on
plot(t, LU*abs(r_lin(1,:)-r(1,:)), 'DisplayName', 'x', 'Color', 'b', 'LineWidth', 1)
plot(t, LU*abs(r_lin(2,:)-r(2,:)), 'DisplayName', 'y', 'Color', 'r', 'LineWidth', 1)
plot(t, LU*abs(r_lin(3,:)-r(3,:)), 'DisplayName', 'z', 'Color', 'g', 'LineWidth', 1)
legend show; grid on;
xlabel('t [TU]'); ylabel('Coordinate Residuals [km]');
fontsize(16,'points'); axis tight

% =========================================================================
%% SECTION 5: QPO Frequency Grid Construction
% =========================================================================
% Generate the set of integer frequency combinations J_QPO that define the
% Fourier basis for the QPO torus. Each row j of J_QPO is an integer pair
% [j1, j2] such that:
%
%   omega_j = j1*w0 + j2*w_P
%
% subject to:
%   0 <= omega_j <= Bh    (bandwidth constraint)
%   |j2| <= max_order_wP  (limits order of quasi-periodic harmonics)
%
% The bandwidth Bh = M*w0 + 2*w_P matches the paper's choice, retaining
% all combinations up to order M in the PO frequency and order 2 in w_P.
%
% Initial Fourier coefficients for the torus are seeded from the PO
% decomposition: QPO-specific frequencies (those with j2 != 0) are
% zeroed out, since the PO itself lives at epsilon = 0.
% =========================================================================

Bh           = M*w_fund(1) + 2*w_P;    % frequency bandwidth [rad/TU]
max_order_wP = 2;                       % maximum order of w_P harmonics

w_fund_QPO = [w_fund(1); w_P];          % two-frequency basis vector
J_QPO      = Freq_Comb(Bh, w_fund_QPO, max_order_wP);
ds         = 5e-4;                      % initial center manifold step size [LU]

% Build the full two-frequency basis
w_fund = [2*pi/Tp, w_P].';
Jx = J_QPO; Jy = J_QPO; Jz = J_QPO;

% Compute frequencies; zero out any combination involving w_P
% (those coefficients are initialized to zero since PO has no w_P content)
freq = Jx * w_fund;
freq(Jx(:,2) ~= 0) = 0;

% Fit initial Fourier coefficients from the PO trajectory
X   = [cos(freq*t); sin(freq*t)];
bin = pinv(X*X.') * X * r(1,:).';
Ax  = bin(1:length(Jx)).';   Bx = bin(length(Jx)+1:end).';

% Consolidate all DC (zero-frequency) components into a single term
Ax(freq==0) = [sum(Ax(freq==0)), zeros(1, length(freq(freq==0))-1)];
Bx(freq==0) = [sum(Bx(freq==0)), zeros(1, length(freq(freq==0))-1)];

Y   = [cos(freq*t); sin(freq*t)];
bin = pinv(Y*Y.') * Y * r(2,:).';
Ay  = bin(1:length(Jy)).';   By = bin(length(Jy)+1:end).';
Ay(freq==0) = [sum(Ay(freq==0)), zeros(1, length(freq(freq==0))-1)];
By(freq==0) = [sum(By(freq==0)), zeros(1, length(freq(freq==0))-1)];

Z   = [cos(freq*t); sin(freq*t)];
bin = pinv(Z*Z.') * Z * r(3,:).';
Az  = bin(1:length(Jz)).';   Bz = bin(length(Jz)+1:end).';
Az(freq==0) = [sum(Az(freq==0)), zeros(1, length(freq(freq==0))-1)];
Bz(freq==0) = [sum(Bz(freq==0)), zeros(1, length(freq(freq==0))-1)];

A     = [{Ax}; {Ay}; {Az}];
B     = [{Bx}; {By}; {Bz}];
x_old = [Ax, Ay, Az, Bx, By, Bz, w_fund.'];   % flat design vector: [A; B; w]

% =========================================================================
%% SECTION 6: Initial QPO Refinement (First QPO at epsilon = ds)
% =========================================================================
% Solve a constrained nonlinear least-squares problem to find the first
% QPO along the center manifold direction at a prescribed step ds from
% the reference PO.
%
% DESIGN VECTOR x = [Ax, Ay, Az, Bx, By, Bz, w0, wP, phi]
%   - Ax/Ay/Az, Bx/By/Bz : Fourier cosine/sine coefficients (length M each)
%   - w0, wP              : fundamental frequencies (fixed here, pm=[0,0])
%   - phi                 : center manifold phase angle
%
% OBJECTIVE (QPO_fun_first):
%   Minimize Jacobi constant variation: C(t) - C(t=0) = 0 over 1000 points.
%   A dynamically consistent torus has constant Jacobi constant everywhere.
%
% CONSTRAINTS (QPO_con_first):
%   (1) Equality: Collocation -- Fourier trajectory matches ode113
%       propagation on each arc segment
%   (2) Equality: Center manifold displacement -- dot(IC_Q - ICp, xi) = ds
%       where xi is the center eigenvector rotated by phase phi
%   (3) Bounds: B_{x,y,z}(freq=0) = 0 (DC sine terms must vanish)
% =========================================================================

init = [x_old, 0];                          % append phase variable phi

% Lower/upper bounds on design vector
lb = -1.5 * ones(1, length(init));
lb([[3,4,5]*length(J_QPO)+1]) = 0;          % enforce B=0 at zero frequency
ub = abs(lb);

% Collocation arc setup
t_LS     = 0 : Tp/2 : 15*Tp;               % arc start times [TU]
n_points = 15;                              % sample points per arc

% Frequency bounds: hold both frequencies fixed for the initial QPO
pm = [0, 0];
lb(end-2:end-1) = w_fund_QPO.' .* (1-pm);
ub(end-2:end-1) = w_fund_QPO.' .* (1+pm);

% Solver options: interior-point with parallel function evaluations
LS_options = optimoptions('lsqnonlin', ...
    'FunctionTolerance',      1e-8, ...
    'StepTolerance',          1e-8, ...
    'OptimalityTolerance',    1e-8, ...
    'ConstraintTolerance',    1e-8, ...
    'MaxIterations',          250,  ...
    'MaxFunctionEvaluations', 10000, ...
    'Algorithm',              'interior-point', ...
    'UseParallel',            true, ...
    'Display',                'off');

x = lsqnonlin(@(x) QPO_fun_first(x, t_LS, length(J_QPO), J_QPO, u), ...
              init, lb, ub, [], [], [], [], ...
              @(x) QPO_con_first(x, t_LS, length(J_QPO), J_QPO, u, ICp, vec_center(:,1).', ds, n_points), ...
              LS_options);

% Unpack solution
A_Q_x = x(1:length(J_QPO));
A_Q_y = x(length(J_QPO)+1       : 2*length(J_QPO));
A_Q_z = x(2*length(J_QPO)+1     : 3*length(J_QPO));
B_Q_x = x(3*length(J_QPO)+1     : 4*length(J_QPO));
B_Q_y = x(4*length(J_QPO)+1     : 5*length(J_QPO));
B_Q_z = x(5*length(J_QPO)+1     : 6*length(J_QPO));
w_Q   = x(end-2:end-1).';
phi   = x(end);

% Store as cell arrays (indexed by continuation step)
J_Q = [{J_QPO}; {J_QPO}; {J_QPO}];
A_Q = [{A_Q_x}; {A_Q_y}; {A_Q_z}];
B_Q = [{B_Q_x}; {B_Q_y}; {B_Q_z}];

% -------------------------------------------------------------------------
% Compute arclength continuation direction for the first step
% -------------------------------------------------------------------------
% The continuation direction xi is the normalized step taken in Fourier
% coefficient space from the PO embedding to the first QPO.
%
% A weighted version (coeff_xi_w) excludes the mean-frequency (DC)
% coefficients and the frequencies themselves from the arclength metric,
% matching the definition of F in Eq. (12) of the paper. This prevents
% the DC terms (which are constrained separately) from dominating the
% step direction.
% -------------------------------------------------------------------------

x = x(1:end-1);                            % remove phi from design vector
dx = (x - x_old);

% Indices of DC coefficients and frequency variables (excluded from F)
indices = [[0:2]*length(J_QPO)+1, length(x)-[1,0]];
ind     = 1:1:size(x,2);

% Build diagonal weight matrix: zero weight on DC and frequency entries
coeff_xi = zeros(size(dx));
weight   = ones(size(ind));
weight(indices) = [0, 0, 0, 0, 0];         % zero out DC + frequency slots
weight          = diag(weight);

dx_edit = (weight * dx(1,:).');
dstep   = norm(dx(1,:));

coeff_xi(1,:)   = 0;
coeff_xi(1,ind) = dx / norm(dx);           % unweighted arclength direction
coeff_xi_w(1,:) = dx_edit / norm(dx_edit); % weighted arclength direction (used in continuation)
x_old           = x;

% =========================================================================
%% SECTION 7: QPO Family Continuation (Arclength in Fourier Space)
% =========================================================================
% Extend the initial QPO into a family using arclength continuation in
% Fourier coefficient space (Sec 3.4 of the paper, Eqs. 12-13).
%
% At each step ii, a new QPO is predicted by stepping along the previous
% arclength direction coeff_xi_w, then corrected by solving a constrained
% least-squares problem.
%
% TERMINATION criteria:
%   - Jacobi constant change |dC| exceeds 0.1 (family has diverged far)
%   - Step size falls below dstep_min (no convergence possible)
%   - Resonance counter exceeds 4 (stuck near a resonance)
%
% STEP SIZE ADAPTATION:
%   - Halved if Fourier energy does not increase (QPO amplitude stagnating)
%   - Reduced by 25% near frequency resonance (w0 ~ res * wP)
%
% STORED OUTPUTS:
%   A_Q, B_Q    : cell arrays of Fourier coefficients, one column per QPO
%   w_Q         : [2 x N_QPO] matrix of frequency pairs
%   IC_pert     : [N_QPO x 6] matrix of initial conditions for each QPO
% =========================================================================

% --- Continuation parameters ---
dC        = 0.001;      % initial Jacobi constant change (dummy, triggers while)
dstep_max = 1e-1;       % maximum allowed arclength step size
dstep_min = 1e-6;       % minimum allowed arclength step size
count     = 0;          % resonance counter

% Indices excluding DC and frequency components (used for Fourier energy)
e_ind = ind;
e_ind(indices) = [];
E_old = norm(x_old(1, e_ind), 2)^2;

ii          = 2;        % continuation index (1 is the initial QPO)
ds          = 1e-3;     % arclength step size [nondimensional Fourier space]
dstep(ii)   = ds;
n_points    = 20;       % collocation points per arc (increased from initial)

while max(abs(dC)) < 0.1

    % Retrieve current Fourier basis and coefficients from previous step
    J_QPO = cell2mat(J_Q(1,1));
    w_Q(:,ii) = w_Q(:,ii-1);
    A_Q_x = cell2mat(A_Q(1,ii-1)); A_Q_y = cell2mat(A_Q(2,ii-1)); A_Q_z = cell2mat(A_Q(3,ii-1));
    B_Q_x = cell2mat(B_Q(1,ii-1)); B_Q_y = cell2mat(B_Q(2,ii-1)); B_Q_z = cell2mat(B_Q(3,ii-1));

    % --- Predictor: step along weighted arclength direction ---
    init = [A_Q_x, A_Q_y, A_Q_z, B_Q_x, B_Q_y, B_Q_z, w_Q(:,ii).'];
    init(ind) = init(ind) + dstep(ii) * coeff_xi_w(ii-1, ind);

    % NOTE: frequency is also extrapolated using the previous step's trend.
    % This is a finite-difference predictor not documented in the paper;
    % it improves convergence but can mislead if the previous step was anomalous.
    init(end) = init(end) + dx(ii-1, end);

    % Bounds: PO frequency held fixed; wP allowed to vary by +/-50%
    lb = -1.5 * ones(1, length(init));
    lb([[3,4,5]*length(J_QPO)+1]) = 0;
    ub = abs(lb);
    pm = [0, 0.50];
    lb(end-1:end) = w_Q(:,ii).' .* (1-pm);
    ub(end-1:end) = w_Q(:,ii).' .* (1+pm);

    % Pin wP for the first continuation step to stabilize convergence
    if ii < 3
        lb(end) = init(end);
        ub(end) = init(end);
    end

    % --- Corrector: constrained least-squares solve ---
    x = lsqnonlin(@(x) QPO_fun(x, t_LS, length(J_QPO), J_QPO, u), ...
                  init, lb, ub, [], [], [], [], ...
                  @(x) QPO_con(x, t_LS, length(J_QPO), J_QPO, u, ds, coeff_xi_w(ii-1,:), x_old(ii-1,:), n_points), ...
                  LS_options);

    % Unpack corrected solution
    A_Q_x = x(1:length(J_QPO));
    A_Q_y = x(length(J_QPO)+1   : 2*length(J_QPO));
    A_Q_z = x(2*length(J_QPO)+1 : 3*length(J_QPO));
    B_Q_x = x(3*length(J_QPO)+1 : 4*length(J_QPO));
    B_Q_y = x(4*length(J_QPO)+1 : 5*length(J_QPO));
    B_Q_z = x(5*length(J_QPO)+1 : 6*length(J_QPO));
    w_Q(:,ii) = x(end-1:end).';

    A_Q(:,ii) = [{A_Q_x; A_Q_y; A_Q_z}];
    B_Q(:,ii) = [{B_Q_x; B_Q_y; B_Q_z}];
    J_Q(:,ii) = [{J_QPO}; {J_QPO}; {J_QPO}];

    % --- Evaluate Jacobi constant variation along the torus ---
    % Sample the Fourier model over a long time window to check that
    % the torus is dynamically consistent (C should be nearly constant).
    Np_c = 2^14;
    t_c  = linspace(0, 200, Np_c);
    Q1   = w_Q(1,ii) * t_c;
    Q2   = w_Q(2,ii) * t_c;
    r_Q  = Q_to_q(A_Q(:,ii), B_Q(:,ii), [Q1; Q2], J_Q(:,ii));
    v_Q  = Q_to_v(A_Q(:,ii), B_Q(:,ii), [Q1; Q2], J_Q(:,ii), w_Q(:,ii));
    C(ii)= JacobiConstant([r_Q(:,1).', v_Q(:,1).'], u);
    dC   = C(ii) - JacobiConstant([r_Q.', v_Q.'], u);

    fprintf(["\nFrequency set after iteration " + num2str(ii) + ...
             ":\n w_0 = " + num2str(w_Q(1,ii)) + ...
             ", w_P = "   + num2str(w_Q(2,ii)) + "\n"])

    % --- Update arclength direction for next step ---
    dx(ii,:)    = (x - x_old(ii-1,:));
    dx_edit     = (weight * dx(ii,:).').';
    coeff_xi(ii,:)   = dx(ii,:) / norm(dx(ii,:));
    coeff_xi_w(ii,:) = dx_edit  / norm(dx_edit);

    % Store initial condition for this QPO (used in harmonic addition)
    IC_pert(ii,:) = [r_Q(:,1).', v_Q(:,1).'];

    % --- Fourier energy monotonicity check ---
    % The Fourier energy norm(F)^2 should increase with QPO amplitude.
    % If it does not, the step was too large or the corrector diverged.
    E_new          = norm(x(e_ind),        2)^2;
    E_old          = norm(x_old(ii-1,e_ind),2)^2;
    x_old(ii,:)    = x;
    E_min(ii+1)    = E_new - E_old;

    % --- Nearest resonance order (used for resonance detection) ---
    res = round(w_Q(1,ii) / w_Q(2,ii));

    if E_min(ii+1) < 0 || abs(E_min(ii+1)) < 1e-7
        % Fourier energy not increasing: step too large, reduce and retry
        ds = ds * 0.5;

    elseif abs(res*w_Q(2,ii) - w_Q(1,ii)) < 1e-3
        % Near a rational resonance w0 ~ res*wP: reduce step and flag
        % NOTE: the threshold 1e-3 [rad/TU] may need tuning for other orbits
        disp('Resonance!')
        count = count + 1;
        ds    = ds * 0.75;

    else
        % Successful step: advance continuation index
        ii = ii + 1;
    end

    dstep(ii) = ds;     % store (possibly updated) step size

    % --- Termination checks ---
    if dstep(ii) <= dstep_min
        disp('Minimum step size reached with no convergence. Terminating alg.')
        break;
    elseif count > 4
        disp('Stuck at resonance. Terminating alg.')
        break;
    end

end

% =========================================================================
%% SECTION 8: Systematic Harmonic Addition
% =========================================================================
% After constructing the QPO family, the initial Fourier model may have
% insufficient frequency content for larger-amplitude orbits. This section
% iterates through selected family members and adds missing harmonics.
%
% ALGORITHM (per selected QPO):
%   1. Compute phase-point drift: launch a grid of 50x50 initial phase
%      combinations, propagate each for 1 day under true CR3BP dynamics,
%      and compare against the Fourier model (Eq. 23-24 of the paper).
%   2. If max relative error < e_tol: no harmonics needed.
%   3. Otherwise: run NAFF on the Jacobi constant residual dC(t) to
%      identify frequencies not yet in the model.
%   4. Use CVX/Gurobi to find integer combinations J_add matching those
%      frequencies (Sec 3.5, integer relation fitting).
%   5. Append J_add to J_QPO with zero initial coefficients and re-solve
%      the collocation problem.
%   6. Repeat until error < e_tol or no new harmonics can be found.
%
% RELATIVE ERROR METRIC:
%   Both the outer (pre-addition) and inner (post-addition) error checks
%   normalize by distance to the Moon, consistent with Eq. (24):
%     rel_e = ||r_Fourier - r_CR3BP|| / ||r_CR3BP - r_Moon||
% =========================================================================

e_tol  = 0.0001;    % relative position error tolerance (dimensionless)
M_NAFF = 30;        % number of NAFF frequencies to extract from dC residual

% Build 50x50 grid of phase-point combinations [phi1; phi2] in [0, 2pi]^2
phi = linspace(0, 2*pi, 50);
[phi1, phi2] = meshgrid(phi, phi);
phi = [phi1(:).'; phi2(:).'];

% Propagation window for phase-point error: 1 day converted to TU
dt  = 1*24*3600 / TU;
t_t = linspace(0, dt, 500);

% Select which QPO family members to process (every 4th from index 4 to 36)
% Modify this array to target specific members of the family.
index = 4:4:36;

J_add_total = 0;    % cumulative count of harmonics added across all orbits

for ii = 1:length(index)

    % Retrieve the frequency basis for this orbit
    if ii == 1
        J_QPO = cell2mat(J_Q(1,1));         % first orbit: use second-order basis
    else
        J_QPO = cell2mat(J_Q2(1,ii-1));     % subsequent: inherit previous additions
    end

    C_QPO   = JacobiConstant(IC_pert(index(ii),:), u);
    dC      = 0.001;                        % dummy initializer
    w_Q2(:,ii) = w_Q(:,index(ii));

    % Pull Fourier coefficients for this QPO from the continuation family
    A_Q_x = cell2mat(A_Q(1,index(ii)));
    A_Q_y = cell2mat(A_Q(2,index(ii)));
    A_Q_z = cell2mat(A_Q(3,index(ii)));
    B_Q_x = cell2mat(B_Q(1,index(ii)));
    B_Q_y = cell2mat(B_Q(2,index(ii)));
    B_Q_z = cell2mat(B_Q(3,index(ii)));

    % Pad coefficient arrays with zeros for any harmonics added by prior orbits
    A_Q_x = [A_Q_x, zeros(1,J_add_total)]; A_Q_y = [A_Q_y, zeros(1,J_add_total)]; A_Q_z = [A_Q_z, zeros(1,J_add_total)];
    B_Q_x = [B_Q_x, zeros(1,J_add_total)]; B_Q_y = [B_Q_y, zeros(1,J_add_total)]; B_Q_z = [B_Q_z, zeros(1,J_add_total)];

    % If harmonics were added for a prior orbit, solve the collocation
    % problem first with the inherited (padded) frequency set before
    % checking whether additional harmonics are needed for this orbit.
    if J_add_total ~= 0

        initial = [A_Q_x, A_Q_y, A_Q_z, B_Q_x, B_Q_y, B_Q_z, w_Q2(:,ii).'];
        lb = -1.5*ones(1,length(initial)); ub = abs(lb);
        pm = [0.000, 0.000];
        lb(end-1:end) = (1-pm) .* w_Q2(:,ii).';
        ub(end-1:end) = (1+pm) .* w_Q2(:,ii).';

        x = lsqnonlin(@(x) QPO_fun(x, t_LS, length(J_QPO), J_QPO, u), ...
                      initial, lb, ub, [], [], [], [], ...
                      @(x) QPO_con_added(x, t_LS, length(J_QPO), J_QPO, u, IC_pert(index(ii),:), n_points), ...
                      LS_options);

        A_Q_x = x(1:length(J_QPO));
        A_Q_y = x(length(J_QPO)+1   : 2*length(J_QPO));
        A_Q_z = x(2*length(J_QPO)+1 : 3*length(J_QPO));
        B_Q_x = x(3*length(J_QPO)+1 : 4*length(J_QPO));
        B_Q_y = x(4*length(J_QPO)+1 : 5*length(J_QPO));
        B_Q_z = x(5*length(J_QPO)+1 : 6*length(J_QPO));
        w_Q2(:,ii) = x(end-1:end).';
    end

    % Store updated orbit
    A_Q2(:,ii) = [{A_Q_x}; {A_Q_y}; {A_Q_z}];
    B_Q2(:,ii) = [{B_Q_x}; {B_Q_y}; {B_Q_z}];
    J_Q2(:,ii) = [{J_QPO}; {J_QPO}; {J_QPO}];

    % Evaluate Jacobi constant residual over long time arc
    Np_c = 2^15;
    t_c  = linspace(0, 400, Np_c);
    Q1   = w_Q2(1,ii) * t_c;
    Q2   = w_Q2(2,ii) * t_c;
    r_Q  = Q_to_q(A_Q2(:,ii), B_Q2(:,ii), [Q1; Q2], J_Q2(:,ii));
    v_Q  = Q_to_v(A_Q2(:,ii), B_Q2(:,ii), [Q1; Q2], J_Q2(:,ii), w_Q2(:,ii));
    dC   = C_QPO - JacobiConstant([r_Q.', v_Q.'], u);
    dC_prev = max(abs(dC));

    % --- Compute phase-point drift (pre-harmonic-addition baseline) ---
    % Normalize relative error by distance to Moon (Eq. 24 of the paper).
    for jj = 1:length(phi)
        Q     = w_Q2(:,ii) * t_t + phi(:,jj);
        r_Qt  = Q_to_q(A_Q2(:,ii), B_Q2(:,ii), Q, J_Q2(:,ii));
        v_Qt  = Q_to_v(A_Q2(:,ii), B_Q2(:,ii), Q, J_Q2(:,ii), w_Q2(:,ii));
        IC_qt = [r_Qt(:,1).', v_Qt(:,1).'];

        [~,sqt]    = ode113(@(t,s) CR3BP(t,s,u), t_t, [IC_qt], options);
        dist_m     = sqt(:,1:3);
        dist_m(:,1)= dist_m(:,1) - (1-u);              % offset to Moon position
        rel_e(jj)  = max(vecnorm(r_Qt - sqt(:,1:3).', 2, 1) ./ vecnorm(dist_m.', 2, 1));
    end

    w_Q_temp    = w_Q2(:,ii);
    rel_e_prev  = max(rel_e);

    if max(rel_e) < e_tol
        fprintf('\n--------------------------------------------- \nNo added harmonics needed\n')
    end

    % =====================================================================
    % Inner loop: add harmonics until error is below tolerance
    % =====================================================================
    while max(rel_e) > e_tol

        % Run NAFF on dC to identify frequencies present in the residual
        [~,~,w_c] = NAFF(dC.', t_c, 2, M_NAFF);

        % Determine integer bounds for new combinations
        max_int_1 = max(J_QPO(:,1), [], 1) + 1;
        max_int_2 = max(J_QPO(:,2), [], 1) + 1;

        % Solve integer relation problem to match NAFF frequencies
        [J_add] = QPO_Add_Harmonics(w_c, w_Q2(:,ii), max_int_2, max_int_1, J_QPO);

        if isempty(J_add)
            % No new combinable frequencies found: store current result and exit
            A_Q2(:,ii) = [{A_Q_x}; {A_Q_y}; {A_Q_z}];
            B_Q2(:,ii) = [{B_Q_x}; {B_Q_y}; {B_Q_z}];
            J_Q2(:,ii) = [{J_QPO}; {J_QPO}; {J_QPO}];
            w_Q2(:,ii) = w_Q_temp;
            fprintf('\n--------------------------------------------- \nContinuation completed until no new harmonics can be added\n')
            break;
        else
            fprintf('\nAdding harmonic combinations:\n');
            fprintf('[%d %d]\n', J_add.');
        end

        % Append new harmonics to the frequency basis with zero initial amplitudes
        J_add_total = J_add_total + size(J_add, 1);
        w_Q2(:,ii)  = w_Q(:,index(ii));

        A_Q_x = [A_Q_x, zeros(1,size(J_add,1))]; A_Q_y = [A_Q_y, zeros(1,size(J_add,1))]; A_Q_z = [A_Q_z, zeros(1,size(J_add,1))];
        B_Q_x = [B_Q_x, zeros(1,size(J_add,1))]; B_Q_y = [B_Q_y, zeros(1,size(J_add,1))]; B_Q_z = [B_Q_z, zeros(1,size(J_add,1))];

        J_QPO   = [J_QPO; J_add];
        initial = [A_Q_x, A_Q_y, A_Q_z, B_Q_x, B_Q_y, B_Q_z, w_Q2(:,ii).'];

        % Increase collocation points if problem size demands it
        if length(initial) > length(t_LS) * n_points
            n_points = n_points + 5;
        end

        lb = -1.5*ones(1,length(initial)); ub = abs(lb);
        lb(end-1:end) = (1-pm) .* w_Q2(:,ii).';
        ub(end-1:end) = (1+pm) .* w_Q2(:,ii).';

        x = lsqnonlin(@(x) QPO_fun(x, t_LS, length(J_QPO), J_QPO, u), ...
                      initial, lb, ub, [], [], [], [], ...
                      @(x) QPO_con_added(x, t_LS, length(J_QPO), J_QPO, u, IC_pert(index(ii),:), n_points), ...
                      LS_options);

        A_Q_x = x(1:length(J_QPO));
        A_Q_y = x(length(J_QPO)+1   : 2*length(J_QPO));
        A_Q_z = x(2*length(J_QPO)+1 : 3*length(J_QPO));
        B_Q_x = x(3*length(J_QPO)+1 : 4*length(J_QPO));
        B_Q_y = x(4*length(J_QPO)+1 : 5*length(J_QPO));
        B_Q_z = x(5*length(J_QPO)+1 : 6*length(J_QPO));
        w_Q_temp = x(end-1:end).';

        % Update dC with new Fourier model
        Q1 = w_Q_temp(1)*t_c; Q2 = w_Q_temp(2)*t_c;
        A_Q_temp = [{A_Q_x}; {A_Q_y}; {A_Q_z}];
        B_Q_temp = [{B_Q_x}; {B_Q_y}; {B_Q_z}];
        J_Q_temp = [{J_QPO}; {J_QPO}; {J_QPO}];
        r_Q = Q_to_q(A_Q_temp, B_Q_temp, [Q1; Q2], J_Q_temp);
        v_Q = Q_to_v(A_Q_temp, B_Q_temp, [Q1; Q2], J_Q_temp, w_Q_temp);
        dC  = C_QPO - JacobiConstant([r_Q.', v_Q.'], u);

        % --- Re-evaluate phase-point drift with updated model ---
        % FIX: normalize by distance to Moon (not barycenter distance),
        % consistent with Eq. (24) of the paper and the pre-addition check above.
        for jj = 1:length(phi)
            Q     = w_Q_temp * t_t + phi(:,jj);
            r_Qt  = Q_to_q(A_Q_temp, B_Q_temp, Q, J_Q_temp);
            v_Qt  = Q_to_v(A_Q_temp, B_Q_temp, Q, J_Q_temp, w_Q_temp);
            IC_qt = [r_Qt(:,1).', v_Qt(:,1).'];

            [~,sqt]     = ode113(@(t,s) CR3BP(t,s,u), t_t, [IC_qt], options);
            dist_m      = sqt(:,1:3);
            dist_m(:,1) = dist_m(:,1) - (1-u);          % FIX: offset to Moon position
            rel_e(jj)   = max(vecnorm(r_Qt - sqt(:,1:3).', 2, 1) ./ vecnorm(dist_m.', 2, 1));
        end

        if max(rel_e) < e_tol
            fprintf('\n--------------------------------------------- \nRelative error is now below tolerance\n')
        end

        rel_e_prev = max(rel_e);

    end % while max(rel_e) > e_tol

    % Store final result for this orbit
    A_Q2(:,ii) = [{A_Q_x}; {A_Q_y}; {A_Q_z}];
    B_Q2(:,ii) = [{B_Q_x}; {B_Q_y}; {B_Q_z}];
    J_Q2(:,ii) = [{J_QPO}; {J_QPO}; {J_QPO}];

end % for ii = 1:length(index)

% =========================================================================
%% SECTION 9: Figures
% =========================================================================

% --- Figure: QPO family with second-order harmonics only ---
t_Q = linspace(0, 1000, 10000);

figure()
for ii = 1:length(index)
    scale = LU / 10^5;
    Q     = w_Q(:,index(ii)) * t_Q;
    r_Qt  = Q_to_q(A_Q(:,index(ii)), B_Q(:,index(ii)), Q, J_Q(:,index(ii))) * scale;

    nexttile
    hold on
    line(r_Qt(1,:), r_Qt(2,:), r_Qt(3,:), 'Color', [1, 0, 0, 0.05], 'LineWidth', 2)
    plot3(scale*r(1,:), scale*r(2,:), scale*r(3,:), 'b', 'LineWidth', 2)
    opts.Position = scale * [1-u; 0; 0];
    opts.Units    = 'EM_km_scale';
    planet3D('Moon', opts);
    xlabel('x [\times10^5 km]'); ylabel('y [\times10^5 km]'); zlabel('z [\times10^5 km]');
    fontsize(16,'points'); axis equal; axis tight; grid on; view(-50, 15)
end

% --- Figure: QPO family with added harmonics ---
t_Q = linspace(0, 1000, 10000);
figure()
for ii = 1:length(index)
    scale = LU / 10^5;
    Q     = w_Q2(:,ii) * t_Q;
    r_Qt  = Q_to_q(A_Q2(:,ii), B_Q2(:,ii), Q, J_Q2(:,ii)) * scale;

    nexttile
    hold on
    line(r_Qt(1,:), r_Qt(2,:), r_Qt(3,:), 'Color', [1, 0, 0, 0.05], 'LineWidth', 2)
    plot3(scale*r(1,:), scale*r(2,:), scale*r(3,:), 'b', 'LineWidth', 2)
    opts.Position = scale * [1-u; 0; 0];
    opts.Units    = 'EM_km_scale';
    planet3D('Moon', opts);
    xlabel('x [\times10^5 km]'); ylabel('y [\times10^5 km]'); zlabel('z [\times10^5 km]');
    fontsize(16,'points'); axis equal; axis tight; grid on; view(-50, 15)
end

% =========================================================================
%% SECTION 10: Phase-Point Error Calculation (Final Validation)
% =========================================================================
% Compute one-day drift errors for both the 2nd-order and added-harmonic
% models over a 50x50 grid of initial phase combinations.
%
% Four error metrics are recorded per phase-point (Eqs. 23-24):
%   e_pos  : absolute position error [km]
%   rel_e  : relative position error, normalized by distance to Moon [-]
%   e_vel  : absolute velocity error [LU/TU]
%   rel_v  : relative velocity error, normalized by velocity magnitude [-]
% =========================================================================

phi = linspace(0, 2*pi, 50);
[phi1, phi2] = meshgrid(phi, phi);
phi = [phi1(:).'; phi2(:).'];

dt  = 1*24*3600 / TU;
t_t = linspace(0, dt, 500);

for ii = 1:length(phi)
    for jj = 1:length(index)

        % --- Added-harmonics model ---
        Q    = w_Q2(:,jj) * t_t + phi(:,ii);
        r_Qt = Q_to_q(A_Q2(:,jj), B_Q2(:,jj), Q, J_Q2(:,jj));
        v_Qt = Q_to_v(A_Q2(:,jj), B_Q2(:,jj), Q, J_Q2(:,jj), w_Q2(:,jj));
        IC_qt= [r_Qt(:,1).', v_Qt(:,1).'];

        [~,sqt]    = ode113(@(t,s) CR3BP(t,s,u), t_t, [IC_qt], options);
        r_moon     = sqt(:,1:3);
        r_moon(:,1)= r_moon(:,1) - (1-u);      % position relative to Moon

        e_pos2(jj,ii) = LU * max(vecnorm(r_Qt - sqt(:,1:3).', 2, 1), [], 2);
        rel_e2(jj,ii) = max(vecnorm(r_Qt - sqt(:,1:3).', 2, 1) ./ vecnorm(r_moon.', 2, 1));
        rel_v2(jj,ii) = max(vecnorm(v_Qt - sqt(:,4:6).', 2, 1) ./ vecnorm(sqt(:,4:6).', 2, 1));
        e_vel2(jj,ii) = max(vecnorm(v_Qt - sqt(:,4:6).', 2, 1), [], 2);

        % --- Second-order harmonics model (baseline comparison) ---
        Q    = w_Q(:,index(jj)) * t_t + phi(:,ii);
        r_Qt = Q_to_q(A_Q(:,index(jj)), B_Q(:,index(jj)), Q, J_Q(:,index(jj)));
        v_Qt = Q_to_v(A_Q(:,index(jj)), B_Q(:,index(jj)), Q, J_Q(:,index(jj)), w_Q(:,index(jj)));
        IC_qt= [r_Qt(:,1).', v_Qt(:,1).'];

        [~,sqt]    = ode113(@(t,s) CR3BP(t,s,u), t_t, [IC_qt], options);
        r_moon     = sqt(:,1:3);
        r_moon(:,1)= r_moon(:,1) - (1-u);

        e_pos(jj,ii) = LU * max(vecnorm(r_Qt - sqt(:,1:3).', 2, 1), [], 2);
        e_vel(jj,ii) = max(vecnorm(v_Qt - sqt(:,4:6).', 2, 1), [], 2);
        rel_e(jj,ii) = max(vecnorm(r_Qt - sqt(:,1:3).', 2, 1) ./ vecnorm(r_moon.', 2, 1));
        rel_v(jj,ii) = max(vecnorm(v_Qt - sqt(:,4:6).', 2, 1) ./ vecnorm(sqt(:,4:6).', 2, 1));

    end
end

% --- Violin plot: relative position error ---
figure()
hold on
vp1 = violinplot(rel_e.',  'FaceColor', 'b', 'DisplayName', '2-Harmonics',    'DensityDirection', 'negative');
vp2 = violinplot(rel_e2.', 'FaceColor', 'r', 'DisplayName', 'Added Harmonics','DensityDirection', 'positive');
plot(1:1:length(index), max(rel_e,  [], 2), 'b .', 'MarkerSize', 10)
plot(1:1:length(index), max(rel_e2, [], 2), 'r .', 'MarkerSize', 10)
xticklabels([string(round(w_Q2(2,:), 4))])
xlabel('\omega_P'); ylabel('(E_\infty)_{rel}');
grid on; axis tight;
legend([vp1(1), vp2(1)], {'2-Harmonics', 'Added Harmonics'})
fontsize(16,'points');
ylim([0, max([rel_e; rel_e2], [], 'all')]);
ax = gca; ax.YAxis.Exponent = 0;

% --- Violin plot: relative velocity error ---
figure()
hold on
vp1 = violinplot(rel_v.',  'FaceColor', 'b', 'DisplayName', '2-Harmonics',    'DensityDirection', 'negative');
vp2 = violinplot(rel_v2.', 'FaceColor', 'r', 'DisplayName', 'Added Harmonics','DensityDirection', 'positive');
plot(1:1:length(index), max(rel_v,  [], 2), 'b .', 'MarkerSize', 10)
plot(1:1:length(index), max(rel_v2, [], 2), 'r .', 'MarkerSize', 10)
xticklabels([string(round(w_Q2(2,:), 4))])
xlabel('\omega_P'); ylabel('(V_\infty)_{rel}');
grid on; axis tight;
legend([vp1(1), vp2(1)], {'2-Harmonics', 'Added Harmonics'})
fontsize(16,'points');
ylim([0, max([rel_v; rel_v2], [], 'all')]);
ax = gca; ax.YAxis.Exponent = 0;

