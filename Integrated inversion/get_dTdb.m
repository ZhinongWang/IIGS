function dT_db = get_dTdb(fq,vs,vp,rho,thk,c)
% compute the partial derivative of the matrix T with respect to S-wave
% velocity
% freq,c -- frequency and phase velocity
% vs,vp,rho,thk -- parameters of the model
% 2019-4-20,ok

w = 2*pi*fq; 
k = w/c;

ds0_db = 2*c^2/vs^3;

r0 = 1-c*c/(vp*vp);
s0 = 1-c*c/(vs*vs);

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

t = 2 - c^2/(vs^2);
mu = rho*vs^2;
gm = vs.^2/c^2;

dgm_db = 2*vs/c^2;
dt_db = 2*c^2/vs^3;
dmu_db = 2*rho*vs;


if c <= vp
    Ca = cosh(k*r*thk);
    Sa = sinh(k*r*thk);
else
    Ca = cos(k*sqrt(-r0)*thk);
    Sa = 1i*sin(k*sqrt(-r0)*thk);
end
if c<= vs
    Cb = cosh(k*s*thk);
    Sb = sinh(k*s*thk);

    dCb_db = sinh(k*s*thk)*k*thk*ds_db;
    dSb_db = cosh(k*s*thk)*k*thk*ds_db;
else
    Cb = cos(k*sqrt(-s0)*thk);
    Sb = 1i*sin(k*sqrt(-s0)*thk);

    dCb_db = 0.5*sin(k*sqrt(-s0)*thk)*k*thk*ds0_db/sqrt(-s0);
    dSb_db = -1i*0.5*cos(k*sqrt(-s0)*thk)*k*thk*ds0_db/sqrt(-s0);
end

dCaCb_db = Ca*dCb_db;
dCaSb_db = Ca*dSb_db;
dSaCb_db = Sa*dCb_db;
dSaSb_db = Sa*dSb_db;

dT_db = zeros(5);
Q = zeros(5,1);
dQ_db = zeros(5,1);

for m = 0:4
    temp = t^m/r/s + 2^m*r*s;
    Q(m+1) = temp*Sa*Sb;
    dQ_db(m+1) = (1/r/s*(m*t^(m-1)*dt_db - t^m*ds_db/s) + 2^m*r*ds_db)*Sa*Sb + temp*Sa*dSb_db;
end

dT_db(1,1) = 2*gm*dgm_db*(-4*t+(t^2+4)*Ca*Cb-Q(3)) +...
             gm^2*(-4*dt_db+2*t*dt_db*Ca*Cb+(t^2+4)*Ca*dCb_db-dQ_db(3));
dT_db(1,2) = ((2*gm*dgm_db*mu-gm^2*dmu_db)/mu^2)*((2+t)*(1-Ca*Cb)+Q(2)) +...
             gm^2/mu*(dt_db*(1-Ca*Cb)-(2+t)*Ca*dCb_db+dQ_db(2));
dT_db(1,3) = 1/mu*((dgm_db-gm*dmu_db/mu)*(Ca*Sb/s-r*Sa*Cb) +...
             gm*(Ca/s*(dSb_db-Sb*ds_db/s)-r*dSaCb_db));
dT_db(1,4) = 1/mu*((dgm_db-gm*dmu_db/mu)*(s*Ca*Sb-Sa*Cb/r) +...
             gm*(Ca*(ds_db*Sb+s*dSb_db)-dSaCb_db/r));
dT_db(1,5) = 1/mu^2*(2*gm*(dgm_db-gm*dmu_db/mu)*(2*(1-Ca*Cb)+Q(1)) + gm^2*(-2*Ca*dCb_db+dQ_db(1)));

dT_db(2,1) = gm*((dmu_db*gm+2*mu*dgm_db)*(-2*t*(t+2)*(1-Ca*Cb)-Q(4)) +...
             mu*gm*(-4*dt_db*(t+1)*(1-Ca*Cb)+2*t*(t+2)*Ca*dCb_db-dQ_db(4)));
dT_db(2,1) = 2*dT_db(2,1);
dT_db(2,2) = dCaCb_db - dT_db(1,1);
dT_db(2,2) = 2*dT_db(2,2);
dT_db(2,3) = dgm_db*(t/s*Ca*Sb-2*r*Sa*Cb) + gm*(dt_db*Ca*Sb/s + t/s^2*(dCaSb_db*s-Ca*Sb*ds_db)-2*r*dSaCb_db);
dT_db(2,3) = 2*dT_db(2,3);
dT_db(2,4) = dgm_db*(2*s*Ca*Sb-t/r*Sa*Cb) +...
             gm*(2*(ds_db*Ca*Sb+s*dCaSb_db) - dt_db*Sa*Cb/r - t/r*dSaCb_db);
dT_db(2,4) = 2*dT_db(2,4);
dT_db(2,5) = 2*dT_db(1,2);

dT_db(3,1) = (dmu_db*gm+mu*dgm_db)*(4*s*Ca*Sb-t^2/r*Sa*Cb) +...
             mu*gm*(4*(ds_db*Ca*Sb+s*dCaSb_db) - 2*t*dt_db/r*Sa*Cb - t^2/r*dSaCb_db);
dT_db(3,2) = -dT_db(2,4)*0.5;
dT_db(3,3) = dCaCb_db;
dT_db(3,4) = -Sa/r*(ds_db*Sb+s*dSb_db);
dT_db(3,5) = -dT_db(1,4);

dT_db(4,1) = (dmu_db*gm+mu*dgm_db)*(t^2/s*Ca*Sb-4*r*Sa*Cb) +...
             mu*gm*(2*t*dt_db*Ca*Sb/s + t^2/s^2*(Ca*dSb_db*s-Ca*Sb*ds_db) - 4*r*Sa*dCb_db);
dT_db(4,2) = -dT_db(2,3)*0.5;
dT_db(4,3) = -r*Sa/s^2*(dSb_db*s - Sb*ds_db);
dT_db(4,4) = dT_db(3,3);
dT_db(4,5) = -dT_db(1,3);

dT_db(5,1) = 2*mu*gm*(dmu_db*gm+mu*dgm_db)*(8*t^2*(1-Ca*Cb)+Q(5)) +...
             mu^2*gm^2*(8*t*(2*dt_db*(1-Ca*Cb)-t*dCaCb_db)+dQ_db(5));
dT_db(5,2) = dT_db(2,1)*0.5;
dT_db(5,3) = -dT_db(4,1);
dT_db(5,4) = -dT_db(3,1);
dT_db(5,5) = dT_db(1,1);