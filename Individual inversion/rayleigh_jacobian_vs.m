function Jac = rayleigh_jacobian_vs(vs, bos, rho, thk)
global MAXROOT THEDATA VRN

Jac = zeros(VRN, length(vs));
for k = 1:MAXROOT
    fqn = length(THEDATA{k}.freq);
    
    dfdc = zeros(fqn,1);
    dfdvs = zeros(fqn,length(vs));
    dcdvs = zeros(fqn, length(vs));
    for fi = 1:fqn
        dfdc(fi) = GetDFDC(THEDATA{k}.freq(fi), THEDATA{k}.vr(fi), vs, bos, rho, thk);
    end
    
    for j = 1:length(vs)
        for fi = 1:length(THEDATA{k}.freq)
            if abs(dfdc(fi)) <= 1e-6, continue;end
            dfdvs(fi,j) = GetDFDVS(j, THEDATA{k}.freq(fi), THEDATA{k}.vr(fi), vs, bos, rho, thk);
            dcdvs(fi,j) = -dfdvs(fi,j)/dfdc(fi);       
        end
        for i = THEDATA{k}.ind(1):THEDATA{k}.ind(end)
            Jac(i,j) = dcdvs(i-THEDATA{k}.ind(1)+1,j);  
        end
    end
end