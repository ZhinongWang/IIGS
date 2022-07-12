function [v,dv_db] = get_dvdb(vs,vp,rho,c)
% input:
% vs,vp,rho - S-wave velocity, P-wave velocity, and density of the
% half-space
% c -- phase velocity
% output:
% v -- the vector V of the half-space
% dv_db -- the partial derivative of vector V with respect to S-wave
% velocity
% 2019-4-20, ok

v = zeros(5,1);
dv_db = zeros(5,1);

r0 = 1 - c^2/vp^2;
s0 = 1 - c^2/vs^2;
ds0_db = 2*c^2/vs^3;

mu = rho*vs^2;
t = 2 - c^2/vs^2;

dmu_db = 2*rho*vs;
dt_db = ds0_db;

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
ds_db = 0.5*ds0_db/s;

v(1) = 1-r*s;
v(2) = 2*mu*(t-2*r*s);
v(3) = mu*s*(2-t);
v(4) = -mu*r*(2-t);
v(5) = mu^2*(4*r*s-t^2);

dv_db(1) = -r*ds_db;
dv_db(2) = 2*(dmu_db*(t-2*r*s) + mu*(dt_db-2*r*ds_db));
dv_db(3) = (dmu_db*s+mu*ds_db)*(2-t) - mu*s*dt_db;
dv_db(4) = -r*(dmu_db*(2-t)-mu*dt_db);
dv_db(5) = mu*(2*dmu_db*(4*r*s-t*t) + mu*(4*r*ds_db-2*t*dt_db));

