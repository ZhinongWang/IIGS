function [rdc,rdl] = rayleigh_forward_vsvp(vs, vp, bos, rho, thk)

global EXPDATA THEDATA MAXROOT MYDATA_THE MYDATA_EXP VRN

guidedphase2(unique(MYDATA_EXP(:,2)), vs, vp, bos, rho, thk);

n = 0;
for k = 1:MAXROOT
    n = n + length(find(THEDATA{k}.vr));
end

s1 = 0; s2 = 0;
for k = 1:VRN
    s2 = s2 + ((MYDATA_THE(k,3) - MYDATA_EXP(k,3))/MYDATA_EXP(k,3))^2;
    if abs(MYDATA_THE(k,3)) <= 1e-6, continue; end
    s1 = s1 + ((log(MYDATA_THE(k,3)) - log(MYDATA_EXP(k,3)))/log(MYDATA_EXP(k,3)))^2;
end

rdc = sqrt(s1 / n);
rdl = sqrt(s2 / VRN);



