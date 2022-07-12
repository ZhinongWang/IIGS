function preprocess
global VRN EXPDATA THEDATA MAXROOT MYDATA_EXP MYDATA_THE Wd 
VRN = length(MYDATA_EXP(:,1));   
MAXROOT = max(MYDATA_EXP(:,1)); 
MYDATA_THE = MYDATA_EXP;
EXPDATA = cell(MAXROOT,1);       
for j = 1:MAXROOT
    ind = find(MYDATA_EXP(:,1) == j);
    EXPDATA{j}.mode = MYDATA_EXP(ind,1);
    EXPDATA{j}.freq = MYDATA_EXP(ind,2);
    EXPDATA{j}.vr = MYDATA_EXP(ind,3);
    EXPDATA{j}.ind = [ind(1) ind(end)]; 
end
THEDATA = EXPDATA;
Wd = diag(1./log(MYDATA_EXP(:,3)));
end





