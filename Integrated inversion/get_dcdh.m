function dcdh = get_dcdh(fq,vr,vs,vp,rho,thk,tag)
% compute the partial derivatives of the phase velocity with respect to
% layer thickness
% freq,vr -- frequency and phase velocity
% vs,vp,rho,thk -- parameters of the model
% tag -- sequence number of the layer

dfdc = reduced_delta_dfdc(fq,vs,vp,rho,thk,vr);
dfdh = reduced_delta_dfdh(fq,vs,vp,rho,thk,vr,tag);
dcdh = -dfdh/dfdc;