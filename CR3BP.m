
function ds = CR3BP(t,S,u)
% CR3BP  Equations of motion for the Circular Restricted Three-Body Problem.
%
%   Implements the standard CR3BP EOM in the rotating synodic frame.
%   The Moon is located at [1-mu, 0, 0] and the Earth at [-mu, 0, 0].
%
%   INPUTS:
%     t   - current time [TU] (unused, autonomous system)
%     S   - [6x1] state vector [x, y, z, xdot, ydot, zdot]
%     u   - CR3BP mass parameter mu
%
%   OUTPUT:
%     ds  - [6x1] time derivative of S

x=S(1); y=S(2); z=S(3); xdot=S(4); ydot=S(5); zdot=S(6);

r1 = sqrt((x+u)^2   + y^2 + z^2);  % distance from Earth to spacecraft
r2 = sqrt((x+u-1)^2 + y^2 + z^2);  % distance from Moon  to spacecraft

ds(1:3) = [xdot, ydot, zdot];
ds(4)   =  2*ydot + x - (1-u)*(x+u)/r1^3   - u*(x+u-1)/r2^3;
ds(5)   = -2*xdot + y - (1-u)*(y)/r1^3     - u*(y)/r2^3;
ds(6)   =             - (1-u)*(z)/r1^3     - u*(z)/r2^3;
ds      = ds';
end
