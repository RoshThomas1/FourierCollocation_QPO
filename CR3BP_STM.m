function ds = CR3BP_STM(t,S,u)
% CR3BP_STM  Equations of motion for the CR3BP with State Transition Matrix.
%
%   Propagates the 6-state trajectory and the 6x6 STM simultaneously,
%   yielding a 42-element state vector [x,y,z,xd,yd,zd, phi(1:36)].
%
%   The STM satisfies: d/dt[Phi] = A(t) * Phi, where A is the
%   time-varying Jacobian of the CR3BP equations of motion.
%
%   INPUTS:
%     t   - current time [TU] (unused, autonomous system)
%     S   - [42x1] state vector: [position; velocity; STM(6x6) reshaped]
%     u   - CR3BP mass parameter mu
%
%   OUTPUT:
%     ds  - [42x1] time derivative of S

x=S(1); y=S(2); z=S(3); xdot=S(4); ydot=S(5); zdot=S(6);

% Jacobian matrix A = df/dx of the CR3BP EOM (auto-generated symbolic form)
A=[0,0,0,1,0,0;0,0,0,0,1,0;0,0,0,0,0,1;(u-1.0)/((u+x)^2+y^2+z^2)^(3/2)-(1.0*u)/(y^2+z^2+(u+x-1.0)^2)^(3/2)-(1.5*(u-1.0)*(u+x)*(2.0*u+2.0*x))/((u+x)^2+y^2+z^2)^(5/2)+(1.5*u*(u+x-1.0)*(2.0*u+2.0*x-2.0))/(y^2+z^2+(u+x-1.0)^2)^(5/2)+1.0,(3.0*u*y*(u+x-1.0))/(y^2+z^2+(u+x-1.0)^2)^(5/2)-(3.0*y*(u-1.0)*(u+x))/((u+x)^2+y^2+z^2)^(5/2),(3.0*u*z*(u+x-1.0))/(y^2+z^2+(u+x-1.0)^2)^(5/2)-(3.0*z*(u-1.0)*(u+x))/((u+x)^2+y^2+z^2)^(5/2),0,2,0;(1.5*u*y*(2.0*u+2.0*x-2.0))/(y^2+z^2+(u+x-1.0)^2)^(5/2)-(1.5*y*(u-1.0)*(2.0*u+2.0*x))/((u+x)^2+y^2+z^2)^(5/2),(u-1.0)/((u+x)^2+y^2+z^2)^(3/2)-(1.0*u)/(y^2+z^2+(u+x-1.0)^2)^(3/2)+(3.0*u*y^2)/(y^2+z^2+(u+x-1.0)^2)^(5/2)-(3.0*y^2*(u-1.0))/((u+x)^2+y^2+z^2)^(5/2)+1.0,(3.0*u*y*z)/(y^2+z^2+(u+x-1.0)^2)^(5/2)-(3.0*y*z*(u-1.0))/((u+x)^2+y^2+z^2)^(5/2),-2,0,0;(1.5*u*z*(2.0*u+2.0*x-2.0))/(y^2+z^2+(u+x-1.0)^2)^(5/2)-(1.5*z*(u-1.0)*(2.0*u+2.0*x))/((u+x)^2+y^2+z^2)^(5/2),(3.0*u*y*z)/(y^2+z^2+(u+x-1.0)^2)^(5/2)-(3.0*y*z*(u-1.0))/((u+x)^2+y^2+z^2)^(5/2),(u-1.0)/((u+x)^2+y^2+z^2)^(3/2)-(1.0*u)/(y^2+z^2+(u+x-1.0)^2)^(3/2)+(3.0*u*z^2)/(y^2+z^2+(u+x-1.0)^2)^(5/2)-(3.0*z^2*(u-1.0))/((u+x)^2+y^2+z^2)^(5/2),0,0,0];

phis   = S(7:42);
phism  = reshape(phis, 6, 6);
phidot = A * phism;

r1 = sqrt((x+u)^2   + y^2 + z^2);
r2 = sqrt((x+u-1)^2 + y^2 + z^2);

ds(1:3) = [xdot, ydot, zdot];
ds(4)   =  2*ydot + x - (1-u)*(x+u)/r1^3   - u*(x+u-1)/r2^3;
ds(5)   = -2*xdot + y - (1-u)*(y)/r1^3     - u*(y)/r2^3;
ds(6)   =             - (1-u)*(z)/r1^3     - u*(z)/r2^3;
ds(7:42)= reshape(phidot, 36, 1);
ds      = ds';
end
