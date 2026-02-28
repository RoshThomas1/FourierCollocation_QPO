

function v = OrbitVelocity(A,B,w,t)
% OrbitVelocity  Evaluate the velocity of a Fourier orbit (single-frequency).
%
%   Time derivative of ConstructOrbit. Not used in the main collocation
%   workflow, but available for post-processing.

v  = zeros(length(A), length(t));
Ax = cell2mat(A(1)); Ay = cell2mat(A(2)); Az = cell2mat(A(3));
Bx = cell2mat(B(1)); By = cell2mat(B(2)); Bz = cell2mat(B(3));
wx = cell2mat(w(1)).'; wy = cell2mat(w(2)).'; wz = cell2mat(w(3)).';

v(1,:) = -(Ax.*wx.')*sin(wx*t) + (Bx.*wx.')*cos(wx*t);
v(2,:) = -(Ay.*wy.')*sin(wy*t) + (By.*wy.')*cos(wy*t);
v(3,:) = -(Az.*wz.')*sin(wz*t) + (Bz.*wz.')*cos(wz*t);
end

