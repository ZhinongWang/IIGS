function dT_dc = get_dTdc(fq,vs,vp,rho,thk,c)
% compute the partial derivative of the matrix T with respect to phase
% velocity
% freq,c -- frequency and phase velocity
% vs,vp,rho,thk -- parameters of the model
% 2019-4-18, ok

w = 2*pi*fq; 
k = w/c;
mu = rho*vs^2;

dT_dc = zeros(5);

r0 = 1 - c^2/vp^2;
s0 = 1 - c^2/vs^2;
dr0_dc = -2*c/vp^2;
ds0_dc = -2*c/vs^2;

if c <= vp
    r = sqrt(r0);
    Ca = cosh(k*r*thk);
    Sa = sinh(k*r*thk);
    
    dr_dc = 0.5*dr0_dc/sqrt(r0);
    dCa_dc = sinh(k*r*thk)*k*thk*(dr_dc-r/c);
    dSa_dc = cosh(k*r*thk)*k*thk*(dr_dc-r/c);
else
    r = 1i*sqrt(-r0);
    Ca = cos(k*sqrt(-r0)*thk);
    Sa = 1i*sin(k*sqrt(-r0)*thk);
    
    dr_dc = -1i*0.5*dr0_dc/sqrt(-r0);
    dCa_dc = 1i*sin(-1i*k*r*thk)*k*thk*(dr_dc-r/c);
    %dCa_dc = sin(k*sqrt(-r0)*thk)*k*thk*(0.5*dr0_dc/sqrt(-r0) + sqrt(-r0)/c);
    dSa_dc = cos(-1i*k*r*thk)*k*thk*(dr_dc-r/c);
    %dSa_dc = -1i*cos(k*sqrt(-r0)*thk)*k*thk*(0.5*dr0_dc/sqrt(-r0) + sqrt(-r0)/c);
end

if c <= vs
    s = sqrt(s0);
    Cb = cosh(k*s*thk);
    Sb = sinh(k*s*thk);
    
    ds_dc = 0.5*ds0_dc/sqrt(s0);
    dCb_dc = sinh(k*s*thk)*k*thk*(ds_dc-s/c);
    dSb_dc = cosh(k*s*thk)*k*thk*(ds_dc-s/c);
else
    s = 1i*sqrt(-s0);
    Cb = cos(k*sqrt(-s0)*thk);
    Sb = 1i*sin(k*sqrt(-s0)*thk);
    
    ds_dc = -1i*0.5*ds0_dc/sqrt(-s0);
    dCb_dc = 1i*sin(-1i*k*s*thk)*k*thk*(ds_dc-s/c);
    %dCb_dc = sin(k*sqrt(-s0)*thk)*k*thk*(0.5*ds0_dc/sqrt(-s0) + sqrt(-s0)/c);
    dSb_dc = cos(-1i*k*s*thk)*k*thk*(ds_dc-s/c);
    %dSb_dc = -1i*cos(k*sqrt(-s0)*thk)*k*thk*(0.5*ds0_dc/sqrt(-s0) + sqrt(-s0)/c);
end

dCaCb_dc = dCa_dc*Cb+Ca*dCb_dc;
dCaSb_dc = dCa_dc*Sb+Ca*dSb_dc;
dSaCb_dc = dSa_dc*Cb+Sa*dCb_dc;
dSaSb_dc = dSa_dc*Sb+Sa*dSb_dc;

gm = vs^2/c^2;
t = 2 - 1/gm;

dgm_dc = -2*vs^2/c^3;
dt_dc = -2*c/vs^2;

Q = zeros(5,1);
dQ_dc = zeros(5,1);
drs_dc = dr_dc*s + r*ds_dc;
for m = 0:4
    Q(m+1) = (t^m/r/s + 2^m*r*s)*Sa*Sb;
    dQ_dc(m+1) = ((m*t^(m-1)*dt_dc*r*s-t^m*drs_dc)/(r*s)^2 + 2^m*drs_dc)*Sa*Sb + (t^m/r/s + 2^m*r*s)*dSaSb_dc;
end

dT_dc(1,1) = 2*gm*dgm_dc*(-4*t+(t^2+4)*Ca*Cb-Q(3)) + gm^2*(-4*dt_dc+2*t*dt_dc*Ca*Cb+(t^2+4)*dCaCb_dc-dQ_dc(3));
dT_dc(1,2) = 2*gm*dgm_dc/mu*((2+t)*(1-Ca*Cb)+Q(2)) + gm^2/mu*(dt_dc*(1-Ca*Cb)-(2+t)*dCaCb_dc+dQ_dc(2));
dT_dc(1,3) = dgm_dc/mu*(Ca*Sb/s-r*Sa*Cb) + gm/mu*(dCaSb_dc/s-ds_dc/s^2*Ca*Sb-dr_dc*Sa*Cb-r*dSaCb_dc);
dT_dc(1,4) = dgm_dc/mu*(s*Ca*Sb-Sa*Cb/r) + gm/mu*(ds_dc*Ca*Sb+s*dCaSb_dc-dSaCb_dc/r+dr_dc*Sa*Cb/r^2);
dT_dc(1,5) = 2*gm*dgm_dc/mu^2*(2*(1-Ca*Cb)+Q(1)) + gm^2/mu^2*(-2*dCaCb_dc+dQ_dc(1));

dT21_dc = 2*mu*gm*dgm_dc*(-2*t*(t+2)*(1-Ca*Cb)-Q(4)) + mu*gm^2*(-4*(t+1)*dt_dc*(1-Ca*Cb)+2*t*(t+2)*dCaCb_dc-dQ_dc(4));
dT_dc(2,1) = 2*dT21_dc;
dT22_dc = dCaCb_dc - dT_dc(1,1);
dT_dc(2,2) = 2*dT22_dc;
dT23_dc = dgm_dc*(t/s*Ca*Sb-2*r*Sa*Cb) + gm*(dt_dc*Ca*Sb/s + t*(dCaSb_dc*s-Ca*Sb*ds_dc)/s^2 - 2*(dr_dc*Sa*Cb+r*dSaCb_dc));
dT_dc(2,3) = 2*dT23_dc;
dT24_dc = dgm_dc*(2*s*Ca*Sb-t/r*Sa*Cb) + gm*(2*(ds_dc*Ca*Sb + s*dCaSb_dc) - dt_dc*Sa*Cb/r - t*(dSaCb_dc*r-Sa*Cb*dr_dc)/r^2);
dT_dc(2,4) = 2*dT24_dc;
dT26_dc = dT_dc(1,2);
dT_dc(2,5) = 2*dT26_dc;

dT_dc(3,1) = mu*(dgm_dc*(4*s*Ca*Sb-t^2*Sa*Cb/r) + gm*(4*(ds_dc*Ca*Sb+s*dCaSb_dc) - 2*t*dt_dc*Sa*Cb/r - t^2*(dSaCb_dc*r-Sa*Cb*dr_dc)/r^2));
dT_dc(3,2) = -dT24_dc;
dT_dc(3,3) = dCaCb_dc;
dsr_dc = (ds_dc*r - s*dr_dc)/r^2;
dT_dc(3,4) = -(dsr_dc*Sa*Sb + s/r*dSaSb_dc);
dT_dc(3,5) = -dT_dc(1,4);

dT_dc(4,1) = mu*(dgm_dc*(t^2*Ca*Sb/s-4*r*Sa*Cb) + gm*(2*t*dt_dc*Ca*Sb/s + t^2*(dCaSb_dc*s-Ca*Sb*ds_dc)/s^2 - 4*(dr_dc*Sa*Cb+r*dSaCb_dc)));
dT_dc(4,2) = -dT23_dc;
dT_dc(4,3) = -1/s*((dr_dc - r*ds_dc/s)*Sa*Sb + r*dSaSb_dc);
dT_dc(4,4) = dT_dc(3,3);
dT_dc(4,5) = -dT_dc(1,3);

dT_dc(5,1) = 2*mu^2*gm*dgm_dc*(8*t^2*(1-Ca*Cb)+Q(5)) + mu^2*gm^2*(16*t*dt_dc*(1-Ca*Cb)-8*t^2*dCaCb_dc+dQ_dc(5));
dT_dc(5,2) = dT21_dc;
dT_dc(5,3) = -dT_dc(4,1);
dT_dc(5,4) = -dT_dc(3,1);
dT_dc(5,5) = dT_dc(1,1);