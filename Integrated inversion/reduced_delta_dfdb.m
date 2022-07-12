function dfdb = reduced_delta_dfdb(fq,vs,vp,rho,thk,c,tag)
% compute partial derivative of dispersion function F with respect to the S-wave velocity 
% vp -- P-wave velocity; vs -- S-wave velocity; tag -- layer sequence
% number
% rho -- density; thk -- layer thickness; fq -- frequency; c -- phase
% velocity
% 2019-4-20, ok

cn = length(vs);

x = [0 0 0 0 1];
if tag == cn
    for ic = 1:cn-1
        T = get_T(fq,vs(ic),vp(ic),rho(ic),thk(ic),c);
        x = x*T;
    end
    [~,dv_db] = get_dvdb(vs(cn),vp(cn),rho(cn),c);
    dfdb = x*dv_db;
else
    for ic = 1:cn-1
        if ic == tag
            dT_db = get_dTdb(fq,vs(ic),vp(ic),rho(ic),thk(ic),c);
        else
            dT_db = get_T(fq,vs(ic),vp(ic),rho(ic),thk(ic),c);
        end
        x = x*dT_db;
    end
    [v,~] = get_dvdb(vs(cn),vp(cn),rho(cn),c);
    dfdb = x*v;
end








