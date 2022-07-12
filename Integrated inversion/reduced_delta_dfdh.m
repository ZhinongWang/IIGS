function dfdh = reduced_delta_dfdh(fq,vs,vp,rho,thk,c,tag)
% compute partial derivative of dispersion function F with respect to layer thickness 
% vp -- P-wave velocity; vs -- S-wave velocity; tag -- layer sequence
% number
% rho -- density; thk -- layer thickness; fq -- frequency; c -- phase
% velocity
% 2019-4-20, ok

cn = length(vs);

x = [0 0 0 0 1];
for ic = 1:cn-1
    if ic == tag
        dT_dh = get_dTdh(fq,vs(ic),vp(ic),rho(ic),thk(ic),c);
    else
        dT_dh = get_T(fq,vs(ic),vp(ic),rho(ic),thk(ic),c);
    end
    x = x*dT_dh;
end
[v,~] = get_dvdb(vs(cn),vp(cn),rho(cn),c);
dfdh = x*v;









