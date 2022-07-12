function [v,dv_da] = get_dvda(vs,vp,rho,c)
% input:
% vs,vp,rho - S-wave velocity, P-wave velocity, and density of the
% half-space
% c -- phase velocity
% output:
% v -- the vector V of the half-space
% dv_da -- the partial derivative of vector V with respect to P-wave
% velocity
% 2019-4-23, ok

v = zeros(5,1);
dv_da = zeros(5,1);

r0 = 1 - c^2/vp^2;
s0 = 1 - c^2/vs^2;
dr0_da = 2*c^2/vp^3;

mu = rho*vs^2;
t = 2 - c^2/vs^2;

if c <= vp
    r = sqrt(r0);
else
    r = 1i*sqrt(-r0);
end
if c <= vs
    s = sqrt(s0);
else
    s = 1i*sqrt(-s0);
end
dr_da = 0.5*dr0_da/r;

v(1) = 1-r*s;
v(2) = 2*mu*(t-2*r*s);
v(3) = mu*s*(2-t);
v(4) = -mu*r*(2-t);
v(5) = mu^2*(4*r*s-t^2);

dv_da(1) = -dr_da*s;
dv_da(2) = -4*mu*dr_da*s;
dv_da(3) = 0;
dv_da(4) = -mu*dr_da*(2-t);
dv_da(5) = mu^2*4*dr_da*s;