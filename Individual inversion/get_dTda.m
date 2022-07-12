function dT_da = get_dTda(fq,vs,vp,rho,thk,c)
% compute the partial derivative of the matrix T with respect to P-wave
% velocity
% freq,c -- frequency and phase velocity
% vs,vp,rho,thk -- parameters of the model
% 2019-4-23,ok

w = 2*pi*fq; 
k = w/c;

r0 = 1-c*c/(vp*vp);
s0 = 1-c*c/(vs*vs);

dr0_da = 2*c^2/vp^3;

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

t = 2 - c^2/(vs^2);
mu = rho*vs^2;
gm = vs.^2/c^2;

if c <= vp
    Ca = cosh(k*r*thk);
    Sa = sinh(k*r*thk);
    
    dCa_da = sinh(k*r*thk)*k*thk*dr_da;
    dSa_da = cosh(k*r*thk)*k*thk*dr_da;
else
    Ca = cos(k*sqrt(-r0)*thk);
    Sa = 1i*sin(k*sqrt(-r0)*thk);
    
    dCa_da = 0.5*sin(k*sqrt(-r0)*thk)*k*thk*dr0_da/sqrt(-r0);
    dSa_da = -1i*0.5*cos(k*sqrt(-r0)*thk)*k*thk*dr0_da/sqrt(-r0);
end
if c<= vs
    Cb = cosh(k*s*thk);
    Sb = sinh(k*s*thk);
else
    Cb = cos(k*sqrt(-s0)*thk);
    Sb = 1i*sin(k*sqrt(-s0)*thk);
end

dT_da = zeros(5);
Q = zeros(5,1);
dQ_da = zeros(5,1);

for m = 0:4
    temp = t^m/r/s + 2^m*r*s;
    Q(m+1) = temp*Sa*Sb;
    dQ_da(m+1) = (-dr_da*t^m/r^2/s + 2^m*dr_da*s)*Sa*Sb + temp*dSa_da*Sb;
end

dT_da(1,1) = gm^2*((t^2+4)*dCa_da*Cb - dQ_da(3));
dT_da(1,2) = gm^2/mu*(-(2+t)*dCa_da*Cb + dQ_da(2));
dT_da(1,3) = gm/mu*(dCa_da*Sb/s - Cb*(dr_da*Sa + r*dSa_da));
dT_da(1,4) = gm/mu*(s*dCa_da*Sb - (dSa_da*r-Sa*dr_da)*Cb/r^2);
dT_da(1,5) = gm^2/mu^2*(-2*dCa_da*Cb + dQ_da(1));

dT_da(2,1) = mu*gm^2*(2*t*(t+2)*dCa_da*Cb - dQ_da(4));
dT_da(2,1) = 2*dT_da(2,1);
dT_da(2,2) = dCa_da*Cb - dT_da(1,1);
dT_da(2,2) = 2*dT_da(2,2);
dT_da(2,3) = gm*(t/s*dCa_da*Sb - 2*Cb*(dr_da*Sa + r*dSa_da));
dT_da(2,3) = 2*dT_da(2,3);
dT_da(2,4) = gm*(2*s*dCa_da*Sb - t*Cb*(dSa_da*r - Sa*dr_da)/r^2);
dT_da(2,4) = 2*dT_da(2,4);
dT_da(2,5) = 2*dT_da(1,2);

dT_da(3,1) = mu*gm*(4*s*dCa_da*Sb - t^2*Cb*(dSa_da*r - Sa*dr_da)/r^2);
dT_da(3,2) = -dT_da(2,4)*0.5;
dT_da(3,3) = dCa_da*Cb;
dT_da(3,4) = -s*Sb*(dSa_da*r - Sa*dr_da)/r^2;
dT_da(3,5) = -dT_da(1,4);

dT_da(4,1) = mu*gm*(t^2/s*dCa_da*Sb - 4*Cb*(dr_da*Sa + r*dSa_da));
dT_da(4,2) = -dT_da(2,3)*0.5;
dT_da(4,3) = -Sb/s*(dr_da*Sa + r*dSa_da);
dT_da(4,4) = dT_da(3,3);
dT_da(4,5) = -dT_da(1,3);

dT_da(5,1) = mu^2*gm^2*(-8*t^2*dCa_da*Cb + dQ_da(5));
dT_da(5,2) = dT_da(2,1)*0.5;
dT_da(5,3) = -dT_da(4,1);
dT_da(5,4) = -dT_da(3,1);
dT_da(5,5) = dT_da(1,1);