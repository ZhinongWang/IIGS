function Jac = guided_jacobian_vs(m, bos, rho, init_thk)
global MAXROOT THEDATA VRN init_vp

nlayer = length(rho);
vs = exp(real(m(1:nlayer)));
vp = init_vp;
thk = init_thk;
Jac = zeros(VRN, length(m));

for k = 1:MAXROOT
    fqn = length(THEDATA{k}.freq);    
    dcdvs = zeros(fqn, nlayer);   
    %%
    for i = 1:fqn
        for j = 1:nlayer
            dcdvs(i,j) = reduced_delta_get_dcdb(THEDATA{k}.freq(i), ...
                    THEDATA{k}.vr(i), vs, vp, rho, thk, j);  
        end
    end
    %%
    for j = 1:nlayer
        for i = THEDATA{k}.ind(1):THEDATA{k}.ind(end)
            Jac(i,j) = dcdvs(i-THEDATA{k}.ind(1)+1,j);          
        end 
    end
end


