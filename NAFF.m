
function [A,B,w] = NAFF(y,t,window_order,M)
% NAFF  Numerical Analysis of Fundamental Frequencies algorithm.
%
%   Extracts up to M frequencies and their amplitudes from a time series y
%   using the iterative windowed Fourier approach of Laskar (1993).
%
%   ALGORITHM:
%     0. Extract zero-frequency (mean) component via linear least squares.
%     1. Apply Hanning-type window of order window_order to the signal.
%     2. Find FFT peak to initialize frequency guess.
%     3. Refine [w, A, B] simultaneously via Newton iteration on the
%        3x3 system F = [du, delta_Dc, delta_Ds] = 0 (Eq. 32-33).
%     4. Remove extracted component from signal and repeat.
%
%   INPUTS:
%     y            - [1 x N] time series to decompose
%     t            - [1 x N] time vector
%     window_order - order of the Hanning window (typically 1 or 2)
%     M            - maximum number of frequencies to extract
%
%   OUTPUT:
%     A, B - cosine and sine amplitudes [1 x M]
%     w    - frequencies [1 x M] [rad/TU]

Np = length(t); T = t(end);
fs  = Np / T;
w_s = 2*pi * 0.5 * fs * linspace(-1, 1, Np+1);    % FFT frequency axis

% Hanning-type window: chi(t) = c * (1 - cos(2*pi*t/T))^p
window = 2^window_order * factorial(window_order)^2 / factorial(2*window_order) ...
         * (1 - cos(2*pi*t/T)).^window_order;

% --- Step 0: extract zero-frequency component via least-squares ---
A_iter(1) = sum(y) / Np; B_iter(1) = 0; w_iter(1) = 0;
y_win(1,:) = y .* window;
error = 1; tol = 1e-18;

while error > tol
    [Dc,Ds,dDc,dDs,d2Dc,d2Ds] = ShootingParams(y_win(1,:), t, w_iter(1), Np);
    dydA  = cos(w_iter(1)*t);
    gradY = [dydA].';
    y_NLS = A_iter(1);
    df    = pinv(gradY.'*gradY) * (gradY.') * (y(1,:) - y_NLS).';
    A_iter = A_iter + df;
    if norm(df) < 1e-10; break; end
end

% Remove zero-frequency and prepare for higher-frequency extraction
y(2,:)      = ReduceSignal(y(1,:), t, A_iter(1), B_iter(1), w_iter(1));
y_win(2,:)  = y(2,:) .* window;
y_recon     = ReconSig(A_iter, B_iter, w_iter, t);
ii = ii + 1;

% --- Steps 1-M: iterative frequency extraction ---
while 1
    % Initial guess from FFT peak
    [w_iter(ii), A_iter(ii), B_iter(ii)] = ExtractFreq(y_win(ii,:), t, w_s, Np);

    error = 1; k = 0;
    while error > tol
        [Dc,Ds,dDc,dDs,d2Dc,d2Ds] = ShootingParams(y_win(ii,:), t, w_iter(ii), Np);

        % NAFF objective: maximize |Dc + i*Ds| w.r.t. w
        du       = 2/sqrt(Dc^2+Ds^2) * (Dc*dDc + Ds*dDs);
        delta_Dc = 1/Np * (A_iter(ii)*cos(w_iter(ii)*t) + B_iter(ii)*sin(w_iter(ii)*t)) * [cos(w_iter(ii)*t)].' - Dc;
        delta_Ds = 1/Np * (A_iter(ii)*cos(w_iter(ii)*t) + B_iter(ii)*sin(w_iter(ii)*t)) * [sin(w_iter(ii)*t)].' - Ds;

        [d2u,d_delDc,d_delDs,d_delDcdA,d_delDcdB,d_delDsdA,d_delDsdB] = ...
            JacobianComponents(w_iter(ii), A_iter(ii), B_iter(ii), t, Np, Dc, Ds, dDc, dDs, d2Dc, d2Ds);

        F  = [du, delta_Dc, delta_Ds].'; error = norm(F);
        DF = [d2u,       0,           0;
              d_delDc,   d_delDcdA,   d_delDcdB;
              d_delDs,   d_delDsdA,   d_delDsdB];
        df = -pinv(DF) * F;

        w_iter(ii) = w_iter(ii) + df(1);
        A_iter(ii) = A_iter(ii) + df(2);
        B_iter(ii) = B_iter(ii) + df(3);

        if norm(df) < 1e-10 || k > 20; break; end
        k = k + 1;
    end

    if ii > M; break; end

    % Remove this component and continue to next frequency
    y(ii+1,:)   = ReduceSignal(y(ii,:), t, A_iter(ii), B_iter(ii), w_iter(ii));
    y_win(ii+1,:) = y(ii+1,:) .* window;
    y_recon     = ReconSig(A_iter, B_iter, w_iter, t);
    ii = ii + 1;
end

% Clean up any NaN entries from early termination
w_iter(isnan(w_iter)) = [];
A_iter(isnan(w_iter)) = [];
B_iter(isnan(w_iter)) = 0;

A = A_iter; B = B_iter; w = w_iter;
end

