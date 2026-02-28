
function y = Q_to_q(A,B,Q,J)
% Q_to_q  Map from angle coordinates Q to Cartesian position q.
%
%   Evaluates the Fourier torus embedding:
%     q_i(Q) = sum_j [ A_j * cos(J_j . Q) + B_j * sin(J_j . Q) ]
%   for each coordinate i in {x, y, z}.
%
%   INPUTS:
%     A, B - cell arrays {Ax; Ay; Az}, {Bx; By; Bz}
%     Q    - [2 x N] angle vector [Q1; Q2] = [w0*t + phi1; wP*t + phi2]
%     J    - cell array {Jx; Jy; Jz} of integer combination matrices
%
%   OUTPUT:
%     y    - [3 x N] position in rotating frame [LU]

Jx = cell2mat(J(1)); Jy = cell2mat(J(2)); Jz = cell2mat(J(3));
Ax = cell2mat(A(1)); Ay = cell2mat(A(2)); Az = cell2mat(A(3));
Bx = cell2mat(B(1)); By = cell2mat(B(2)); Bz = cell2mat(B(3));

y(1,:) = Ax*cos(Jx*Q) + Bx*sin(Jx*Q);
y(2,:) = Ay*cos(Jy*Q) + By*sin(Jy*Q);
y(3,:) = Az*cos(Jz*Q) + Bz*sin(Jz*Q);
end

