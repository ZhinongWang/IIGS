function Jac = guided_jacobian_vsvp(m, bos, rho, init_thk)
global MYDATA1 MYDATA2 MYDATA EXPDATA1 EXPDATA2 THEDATA1 THEDATA2
global MAXROOT MYDATA_EXP1 MYDATA_EXP2 MYDATA_THE1 MYDATA_THE2 Wd VRN1 VRN2

nlayer = length(rho);
vs = exp(real(m(1:nlayer)));
vp = exp(real(m(nlayer+1:2*nlayer)));
thk = init_thk;
Jac1 = zeros(VRN1, length(m));
Jac2 = zeros(VRN2, length(m));

for k = 1:2
    fqn = length(THEDATA1{k}.freq);    
    dcdvs = zeros(fqn, nlayer);
    dcdvp = zeros(fqn, nlayer);   
    %%
    for i = 1:fqn
        for j = 1:nlayer
            dcdvs(i,j) = reduced_delta_get_dcdb(THEDATA1{k}.freq(i), ...
                    THEDATA1{k}.vr(i), vs, vp, rho, thk, j);  
            dcdvp(i,j) = reduced_delta_get_dcda(THEDATA1{k}.freq(i), ...
                    THEDATA1{k}.vr(i), vs, vp, rho, thk, j);
        end
    end
    %%
    for j = 1:nlayer
        for i = THEDATA1{k}.ind(1):THEDATA1{k}.ind(end)
            Jac1(i,j) = dcdvs(i-THEDATA1{k}.ind(1)+1,j);           
            Jac1(i,nlayer+j) = dcdvp(i-THEDATA1{k}.ind(1)+1,j);  
        end 
    end
end

for k = 1:2
    fqn = length(THEDATA2{k}.freq);    
    dcdvs = zeros(fqn, nlayer);
    dcdvp = zeros(fqn, nlayer);   
    %%
    for i = 1:fqn
        for j = 1:nlayer
            dcdvs(i,j) = reduced_delta_get_dcdb(THEDATA2{k}.freq(i), ...
                    THEDATA2{k}.vr(i), vs, vp, rho, thk, j);  
            dcdvp(i,j) = reduced_delta_get_dcda(THEDATA2{k}.freq(i), ...
                    THEDATA2{k}.vr(i), vs, vp, rho, thk, j);
        end
    end
    %%
    for j = 1:nlayer
        for i = THEDATA2{k}.ind(1):THEDATA2{k}.ind(end)
            Jac2(i,j) = dcdvs(i-THEDATA2{k}.ind(1)+1,j);          
            Jac2(i,nlayer+j) = dcdvp(i-THEDATA2{k}.ind(1)+1,j);  
        end 
    end
end

Jac = [Jac1;Jac2];







