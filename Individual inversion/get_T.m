function T = get_T(fq,vs,vp,rho,thk,c)
% compute the matrix T of any layer
% fq - frequency£»c - phase velocity
% vs - S-wave velocity of the layer£»vp - P-wave velocity of the layer
% rho - density of the layer£»thk - thickness of the layer
% 2019-4-16,ok

w = 2*pi*fq; 
k = w/c;
mu = rho*vs^2;

T = zeros(5);

r0 = 1 - c*c/(vp*vp);
s0 = 1 - c*c/(vs*vs);
if c <= vp
    r = sqrt(r0);
    Ca = cosh(k*r*thk);
    Sa = sinh(k*r*thk);
else
    r = 1i*sqrt(-r0);
    Ca = cos(k*sqrt(-r0)*thk);
    Sa = 1i*sin(k*sqrt(-r0)*thk);
end

if c <= vs
    s = sqrt(s0);
    Cb = cosh(k*s*thk);
    Sb = sinh(k*s*thk);
else
    s = 1i*sqrt(-s0);
    Cb = cos(k*sqrt(-s0)*thk);
    Sb = 1i*sin(k*sqrt(-s0)*thk);
end

gm = vs^2/c^2;
t = 2 - c^2/(vs^2);

Q = zeros(5,1);
for m = 0:4
    Q(m+1) = (t^m/r/s + 2^m*r*s)*Sa*Sb;
end

T(1,1) = gm^2*(-4*t+(t^2+4)*Ca*Cb-Q(3));
T(1,2) = gm^2/mu*((2+t)*(1-Ca*Cb)+Q(2));
T(1,3) = gm/mu*(Ca*Sb/s-r*Sa*Cb);
T(1,4) = gm/mu*(s*Ca*Sb-Sa*Cb/r);
T(1,5) = (gm/mu)^2*(2*(1-Ca*Cb)+Q(1));

T21 = mu*gm^2*(-2*t*(t+2)*(1-Ca*Cb)-Q(4));
T(2,1) = 2*T21;
T22 = 1+Ca*Cb-T(1,1);
T(2,2) = 2*T22-1;
T23 = gm*(t/s*Ca*Sb-2*r*Sa*Cb);
T(2,3) = 2*T23;
T24 = gm*(2*s*Ca*Sb-t*Sa*Cb/r);
T(2,4) = 2*T24;
T26 = T(1,2);
T(2,5) = 2*T26;

T(3,1) = mu*gm*(4*s*Ca*Sb-t*t*Sa*Cb/r);
T(3,2) = -T24;
T(3,3) = Ca*Cb;
T(3,4) = -s/r*Sa*Sb;
T(3,5) = -T(1,4);

T(4,1) = mu*gm*(t*t/s*Ca*Sb-4*r*Sa*Cb);
T(4,2) = -T23;
T(4,3) = -r/s*Sa*Sb;
T(4,4) = T(3,3);
T(4,5) = -T(1,3);

T(5,1) = mu^2*gm^2*(8*t^2*(1-Ca*Cb)+Q(5));
T(5,2) = T21;
T(5,3) = -T(4,1);
T(5,4) = -T(3,1);
T(5,5) = T(1,1);