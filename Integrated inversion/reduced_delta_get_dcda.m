function dcda = reduced_delta_get_dcda(fq,vr,vs,vp,rho,thk,tag)
% compute partial derivative of phase velocity vr with respect to the P-wave velocity 
% vp -- P-wave velocity; vs -- S-wave velocity; tag -- layer sequence
% number
% rho -- density; thk -- layer thickness; fq -- frequency; vr -- phase
% velocity

dfdc = reduced_delta_dfdc(fq,vs,vp,rho,thk,vr);
dfda = reduced_delta_dfda(fq,vs,vp,rho,thk,vr,tag);
dcda = -dfda/dfdc;