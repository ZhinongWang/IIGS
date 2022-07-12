function dcdb = reduced_delta_get_dcdb(fq,vr,vs,vp,rho,thk,tag)
% compute partial derivative of phase velocity vr with respect to the S-wave velocity 
% vp -- P-wave velocity; vs -- S-wave velocity; tag -- layer sequence
% number
% rho -- density; thk -- layer thickness; fq -- frequency; vr -- phase
% velocity

dfdc = reduced_delta_dfdc(fq,vs,vp,rho,thk,vr);
dfdb = reduced_delta_dfdb(fq,vs,vp,rho,thk,vr,tag);
dcdb = -dfdb/dfdc;