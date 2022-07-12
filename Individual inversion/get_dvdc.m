function [v,dv_dc] = get_dvdc(vs,vp,rho,c)
% input:
% vs,vp,rho - S-wave velocity, P-wave velocity, and density of the
% half-space
% c -- phase velocity
% output:
% v -- the vector V of the half-space
% dv_dc -- the partial derivative of vector V with respect to phase
% velocity
% 2019-4-23, ok

v = zeros(5,1);
dv_dc = zeros(5,1);

r0 = 1 - c^2/vp^2;
s0 = 1 - c^2/vs^2;
dr0_dc = -2*c/vp^2;
ds0_dc = -2*c/vs^2;

mu = rho*vs^2;

gm = vs^2/c^2;
t = 2 - 1/gm;
dt_dc = -2*c/vs^2;

if c <= vp
    r = sqrt(r0);
    dr_dc = 0.5*dr0_dc/r;
else
    r = 1i*sqrt(-r0);
    dr_dc = -1i*0.5*dr0_dc/sqrt(-r0);
end
if c <= vs
    s = sqrt(s0);
    ds_dc = 0.5*ds0_dc/s;
else
    s = 1i*sqrt(-s0);
    ds_dc = -1i*0.5*ds0_dc/sqrt(-s0);
end

v(1) = 1-r*s;
v(2) = 2*mu*(t-2*r*s);
v(3) = mu*s*(2-t);
v(4) = -mu*r*(2-t);
v(5) = mu^2*(4*r*s-t^2);

dv_dc(1) = -(dr_dc*s+r*ds_dc);
dv_dc(2) = 2*mu*(dt_dc - 2*(dr_dc*s + r*ds_dc));
dv_dc(3) = mu*(2*ds_dc-ds_dc*t-s*dt_dc);
dv_dc(4) = -mu*(2*dr_dc-dr_dc*t-r*dt_dc);
dv_dc(5) = mu^2*(4*(dr_dc*s+r*ds_dc)-2*t*dt_dc);