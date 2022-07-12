function preprocess
global MYDATA1 MYDATA2 MYDATA EXPDATA1 EXPDATA2 THEDATA1 THEDATA2
global MAXROOT MYDATA_EXP1 MYDATA_EXP2 MYDATA_THE1 MYDATA_THE2 Wd VRN1 VRN2

VRN2 = length(MYDATA_EXP2(:,1));   
MAXROOT = max(MYDATA_EXP2(:,1)); 
MYDATA_THE2 = MYDATA_EXP2;
EXPDATA2 = cell(MAXROOT,1);       
for j = 1:MAXROOT
    ind = find(MYDATA_EXP2(:,1) == j);
    EXPDATA2{j}.mode = MYDATA_EXP2(ind,1);
    EXPDATA2{j}.freq = MYDATA_EXP2(ind,2);
    EXPDATA2{j}.vr = MYDATA_EXP2(ind,3);
    EXPDATA2{j}.ind = [ind(1) ind(end)]; 
end
THEDATA2 = EXPDATA2;
Wd = diag(1./log(MYDATA_EXP2(:,3)));
end





