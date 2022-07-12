function dT_dh = get_dTdh(fq,vs,vp,rho,thk,c)
% compute the partial derivative of the matrix T with respect to layer
% thickness
% freq,c -- frequency and phase velocity
% vs,vp,rho,thk -- parameters of the model
% 2019-4-20,ok

w = 2*pi*fq; 
k = w/c;

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

t = 2 - c^2/(vs^2);
mu = rho*vs^2;
gm = vs^2/c^2;

if c <= vp
    Ca = cosh(k*r*thk);
    Sa = sinh(k*r*thk);
    
    dCa_dh = sinh(k*r*thk)*k*r;
    dSa_dh = cosh(k*r*thk)*k*r;
else
    Ca = cos(k*sqrt(-r0)*thk);
    Sa = 1i*sin(k*sqrt(-r0)*thk);
    
    dCa_dh = 1i*sin(-1i*k*r*thk)*k*r;
    dSa_dh = cos(-1i*k*r*thk)*k*r;
end
if c<= vs
    Cb = cosh(k*s*thk);
    Sb = sinh(k*s*thk);

    dCb_dh = sinh(k*s*thk)*k*s;
    dSb_dh = cosh(k*s*thk)*k*s;
else
    Cb = cos(k*sqrt(-s0)*thk);
    Sb = 1i*sin(k*sqrt(-s0)*thk);

    dCb_dh = 1i*sin(-1i*k*s*thk)*k*s;
    dSb_dh = cos(-1i*k*s*thk)*k*s;
end

dCaCb_dh = dCa_dh*Cb + Ca*dCb_dh;
dCaSb_dh = dCa_dh*Sb + Ca*dSb_dh;
dSaCb_dh = dSa_dh*Cb + Sa*dCb_dh;
dSaSb_dh = dSa_dh*Sb + Sa*dSb_dh;

dT_dh = zeros(5);
Q = zeros(5,1);
dQ_dh = zeros(5,1);

for m = 0:4
    temp = t^m/r/s + 2^m*r*s;
    Q(m+1) = temp*Sa*Sb;
    dQ_dh(m+1) = temp*dSaSb_dh;
end

dT_dh(1,1) = gm^2*((t^2+4)*dCaCb_dh-dQ_dh(3));
dT_dh(1,2) = gm^2/mu*(-(2+t)*dCaCb_dh+dQ_dh(2));
dT_dh(1,3) = gm/mu*(dCaSb_dh/s - r*dSaCb_dh);
dT_dh(1,4) = gm/mu*(s*dCaSb_dh - dSaCb_dh/r);
dT_dh(1,5) = gm^2/mu^2*(-2*dCaCb_dh + dQ_dh(1));

dT_dh(2,1) = mu*gm^2*(2*t*(t+2)*dCaCb_dh - dQ_dh(4));
dT_dh(2,1) = 2*dT_dh(2,1);
dT_dh(2,2) = dCaCb_dh - dT_dh(1,1);
dT_dh(2,2) = 2*dT_dh(2,2);
dT_dh(2,3) = gm*(t*dCaSb_dh/s -2*r*dSaCb_dh);
dT_dh(2,3) = 2*dT_dh(2,3);
dT_dh(2,4) = gm*(2*s*dCaSb_dh - t*dSaCb_dh/r);
dT_dh(2,4) = 2*dT_dh(2,4);
dT_dh(2,5) = 2*dT_dh(1,2);

dT_dh(3,1) = mu*gm*(4*s*dCaSb_dh - t^2*dSaCb_dh/r);
dT_dh(3,2) = -dT_dh(2,4)*0.5;
dT_dh(3,3) = dCaCb_dh;
dT_dh(3,4) = -s/r*dSaSb_dh;
dT_dh(3,5) = -dT_dh(1,4);

dT_dh(4,1) = gm*mu*(t^2/s*dCaSb_dh - 4*r*dSaCb_dh);
dT_dh(4,2) = -dT_dh(2,3)*0.5;
dT_dh(4,3) = -r/s*dSaSb_dh;
dT_dh(4,4) = dT_dh(3,3);
dT_dh(4,5) = -dT_dh(1,3);

dT_dh(5,1) = mu^2*gm^2*(-8*t*t*dCaCb_dh + dQ_dh(5));
dT_dh(5,2) = dT_dh(2,1)*0.5;
dT_dh(5,3) = -dT_dh(4,1);
dT_dh(5,4) = -dT_dh(3,1);
dT_dh(5,5) = dT_dh(1,1);