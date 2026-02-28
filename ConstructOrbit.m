
function r = ConstructOrbit(A,B,w,t)
% ConstructOrbit  Evaluate a Fourier orbit parameterized by a single frequency.
%
%   Used for PO spectral decomposition validation (Section 4).
%   Each coordinate uses its own (identical) frequency array.
%
%   INPUTS:
%     A, B - cell arrays {Ax; Ay; Az}, {Bx; By; Bz} of Fourier amplitudes
%     w    - cell array {wx; wy; wz} of frequency vectors
%     t    - [1 x N] time vector
%
%   OUTPUT:
%     r    - [3 x N] position trajectory

r  = zeros(length(A), length(t));
Ax = cell2mat(A(1)); Ay = cell2mat(A(2)); Az = cell2mat(A(3));
Bx = cell2mat(B(1)); By = cell2mat(B(2)); Bz = cell2mat(B(3));
wx = cell2mat(w(1)).'; wy = cell2mat(w(2)).'; wz = cell2mat(w(3)).';

r(1,:) = Ax*cos(wx*t) + Bx*sin(wx*t);
r(2,:) = Ay*cos(wy*t) + By*sin(wy*t);
r(3,:) = Az*cos(wz*t) + Bz*sin(wz*t);
end
