function dfda = reduced_delta_dfda(fq,vs,vp,rho,thk,c,tag)
% compute partial derivative of dispersion function F with respect to the P-wave velocity 
% vp -- P-wave velocity; vs -- S-wave velocity; tag -- layer sequence
% number
% rho -- density; thk -- layer thickness; fq -- frequency; c -- phase
% velocity
% 2019-4-23,ok 

cn = length(vs);

x = [0 0 0 0 1];
if tag == cn
    for ic = 1:cn-1
        T = get_T(fq,vs(ic),vp(ic),rho(ic),thk(ic),c);
        x = x*T;
    end
    [~,dv_da] = get_dvda(vs(cn),vp(cn),rho(cn),c);
    dfda = x*dv_da;
else
    for ic = 1:cn-1
        if ic == tag
            dT_da = get_dTda(fq,vs(ic),vp(ic),rho(ic),thk(ic),c);
        else
            dT_da = get_T(fq,vs(ic),vp(ic),rho(ic),thk(ic),c);
        end
        x = x*dT_da;
    end
    [v,~] = get_dvda(vs(cn),vp(cn),rho(cn),c);
    dfda = x*v;
end








