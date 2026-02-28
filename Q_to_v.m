
function y = Q_to_v(A,B,Q,J,w)
% Q_to_v  Map from angle coordinates Q to Cartesian velocity v.
%
%   Time derivative of Q_to_q, using the chain rule:
%     v_i = dq_i/dt = sum_j (J_j . w) * [ -A_j * sin(J_j . Q) + B_j * cos(J_j . Q) ]
%
%   INPUTS:
%     A, B - cell arrays of Fourier amplitudes
%     Q    - [2 x N] angle vector
%     J    - cell array of integer combination matrices
%     w    - [2 x 1] fundamental frequency vector [w0; wP]
%
%   OUTPUT:
%     y    - [3 x N] velocity in rotating frame [LU/TU]

Jx = cell2mat(J(1)); Jy = cell2mat(J(2)); Jz = cell2mat(J(3));
Ax = cell2mat(A(1)); Ay = cell2mat(A(2)); Az = cell2mat(A(3));
Bx = cell2mat(B(1)); By = cell2mat(B(2)); Bz = cell2mat(B(3));

y(1,:) = -(Ax.*(Jx*w).')*sin(Jx*Q) + (Bx.*(Jx*w).')*cos(Jx*Q);
y(2,:) = -(Ay.*(Jy*w).')*sin(Jy*Q) + (By.*(Jy*w).')*cos(Jy*Q);
y(3,:) = -(Az.*(Jz*w).')*sin(Jz*Q) + (Bz.*(Jz*w).')*cos(Jz*Q);
end