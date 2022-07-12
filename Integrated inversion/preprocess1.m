function preprocess
global MYDATA1 MYDATA2 MYDATA EXPDATA1 EXPDATA2 THEDATA1 THEDATA2
global MAXROOT MYDATA_EXP1 MYDATA_EXP2 MYDATA_THE1 MYDATA_THE2 Wd VRN1 VRN2

VRN1 = length(MYDATA_EXP1(:,1));   
MAXROOT = max(MYDATA_EXP1(:,1));  
MYDATA_THE1 = MYDATA_EXP1;
EXPDATA1 = cell(MAXROOT,1);       
for j = 1:MAXROOT
    ind = find(MYDATA_EXP1(:,1) == j);
    EXPDATA1{j}.mode = MYDATA_EXP1(ind,1);
    EXPDATA1{j}.freq = MYDATA_EXP1(ind,2);
    EXPDATA1{j}.vr = MYDATA_EXP1(ind,3);
    EXPDATA1{j}.ind = [ind(1) ind(end)]; 
end
THEDATA1 = EXPDATA1;
end





