function [rdc,rdl] = rayleigh_forward_vsvp(vs, vp, bos, rho, thk)
global MYDATA1 MYDATA2 MYDATA EXPDATA1 EXPDATA2 THEDATA1 THEDATA2
global MAXROOT MYDATA_EXP1 MYDATA_EXP2 MYDATA_THE1 MYDATA_THE2 Wd VRN1 VRN2

guidedphase1(unique(MYDATA_EXP1(:,2)), vs, vp, bos, rho, thk);
guidedphase2(unique(MYDATA_EXP2(:,2)), vs, vp, bos, rho, thk);
n = 0;

for k = 1:2
    n = n + length(find(THEDATA1{k}.vr));
    n = n + length(find(THEDATA2{k}.vr));
end

s1 = 0; s2 = 0;
VRN = VRN1+VRN2;

for k = 1:VRN1
    s2 = s2 + ((MYDATA_THE1(k,3) - MYDATA_EXP1(k,3))/MYDATA_EXP1(k,3))^2;
    if abs(MYDATA_THE1(k,3)) <= 1e-6, continue; end
    s1 = s1 + ((log(MYDATA_THE1(k,3)) - log(MYDATA_EXP1(k,3)))/log(MYDATA_EXP1(k,3)))^2;
end

rdc = sqrt(s1 / n);
rdl = sqrt(s2 / VRN1);



