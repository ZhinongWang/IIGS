function [rd] = rayleigh_forward(vs, bos, rho, thk, dv)
global EXPDATA THEDATA
global MAXROOT MYDATA_THE MYDATA_EXP VRN
global Wd  
for k = 1:MAXROOT
    vrt = rayleighphase(EXPDATA{k}.freq, vs, bos, rho, thk, dv, k);
    THEDATA{k}.vr = vrt(:,end);   
    MYDATA_THE(EXPDATA{k}.ind(1):EXPDATA{k}.ind(end),3) = THEDATA{k}.vr;
end
Wd = diag(1./(MYDATA_EXP(:,4)));
s = (MYDATA_EXP(:,3)-MYDATA_THE(:,3))'*(Wd'*Wd)*(MYDATA_EXP(:,3)-MYDATA_THE(:,3));
rd = sqrt(s / VRN); 