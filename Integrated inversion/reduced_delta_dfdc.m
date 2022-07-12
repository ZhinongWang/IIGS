function dfdc = reduced_delta_dfdc(fq,vs,vp,rho,thk,c)
% compute partial derivative of dispersion function F with respect to phase velocity 
% vp -- P-wave velocity; vs -- S-wave velocity;
% rho -- density; thk -- layer thickness; fq -- frequency; c -- phase
% velocity
% 2019-4-17, ok

cn = length(vs);
x_old = [0 0 0 0 1];
y_old = zeros(1,5);

for ic = 1:cn-1
    T = get_T(fq,vs(ic),vp(ic),rho(ic),thk(ic),c);
    dT_dc = get_dTdc(fq,vs(ic),vp(ic),rho(ic),thk(ic),c);
    x_new = x_old*T;
    y_new = y_old*T + x_old*dT_dc;
    x_old = x_new;
    y_old = y_new;
end
[v,dv_dc] = get_dvdc(vs(cn),vp(cn),rho(cn),c);
dfdc = y_new*v + x_new*dv_dc;


