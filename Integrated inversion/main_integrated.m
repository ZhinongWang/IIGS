clc; clear; close all;
global MYDATA1 MYDATA2 MYDATA EXPDATA1 EXPDATA2 THEDATA1 THEDATA2 INIDATA1 INIDATA2
global MAXROOT MYDATA_EXP1 MYDATA_EXP2 MYDATA_THE1 MYDATA_THE2 Wd VRN1 VRN2

%% Read data
R0 = load("R0.mat");
R0 = R0.R0;
R1 = load("R1.mat");
R1 = R1.R1;
G1 = load("G1.mat");
G1 = G1.G1;
G2 = load("G2.mat");
G2 = G2.G2;
MYDATA1 = [R0;R1];
MYDATA2 = [G1;G2];
MYDATA_EXP1 = MYDATA1;
MYDATA_EXP2 = MYDATA2;

%% Pre-processing
preprocess1;
preprocess2;

%% Initial model
true_vs = [270 450 330 540];
true_vp = [900 1300 1100 1800];
true_thk = [5 10 5];

init_vs = [150 250 400 600];
init_vp = [700 900 1500 1900];
init_thk = [5 10 5];
rho = [1.6 1.8 1.6 2.0];
init_v_rat = init_vp./init_vs;
init_bos = (0.5*init_v_rat.^2-1)./(init_v_rat.^2-1);

%% Initial parameter
dv = 1.0;     
rdc = 0.01;    
maxrun = 20;  
init_m = [init_vs init_vp];

%%
tic
init_m = log(init_m(:));
bos = init_bos(:);
rho = rho(:);
nlayer = length(rho);

damp = logspace(-6,2,10);
dm = zeros(length(init_m),length(damp));
vs = zeros(nlayer,length(damp));
vp = zeros(nlayer,length(damp));
rd = zeros(length(damp), 1);

%% 
delta = zeros(length(init_m));
for ri = 1:length(init_m)-1
    delta(ri,ri)= 1.0;
    delta(ri,ri+1)= -1.0;
end

mcur = init_m;
vscur = exp(real(init_m(1:nlayer)));
vpcur = exp(real(init_m(nlayer+1:2*nlayer)));

for iter = 1:10
    v_ratio = vpcur/vscur;
    boscur = (0.5.*v_ratio.^2-1)./(v_ratio.^2-1);
    [~,rdcur] = guided_forward_vsvp(vscur, vpcur, boscur, rho, init_thk);    
    if iter==1
        for k = 1:MAXROOT
            INIDATA1{k}.freq = THEDATA1{k}.freq;
            INIDATA1{k}.vr = THEDATA1{k}.vr;
            INIDATA2{k}.freq = THEDATA2{k}.freq;
            INIDATA2{k}.vr = THEDATA2{k}.vr;
        end
    end
    mpre = mcur;
    rdpre = rdcur;
    
    %% 
    dd1 = log(MYDATA_EXP1(:,3)) - log(MYDATA_THE1(:,3));
    dd2 = log(MYDATA_EXP2(:,3)) - log(MYDATA_THE2(:,3));
    dd = [dd1;dd2];
    dd(isinf(dd)) = 0.;
    Jac = guided_jacobian_vsvp(mcur, boscur, rho, init_thk);       
    Jac(isnan(Jac)) = 0.; 
    Jac(isinf(Jac)) = 0.;
    
    %%
    figure
    imagesc(real(Jac));
    colormap(turbo);
    caxis([0 1.5]);
    colorbar;
    for j=1:2
        line([0,length(mcur)+1],[EXPDATA1{j}.ind(end)+0.5, EXPDATA1{j}.ind(end)+0.5], ...
            'Color','w','linestyle','--','Linewidth',2);
    end
    for j=1:2
        line([0,length(mcur)+1],[EXPDATA1{2}.ind(end)+EXPDATA2{j}.ind(end)+0.5, ...
            EXPDATA1{2}.ind(end)+EXPDATA2{j}.ind(end)+0.5], ...
            'Color','w','linestyle','--','Linewidth',2);
    end
    line([nlayer+0.5,nlayer+0.5],[0,length(dd)+1], ...
        'Color','w','linestyle','--','Linewidth',2);
    line([2*nlayer+0.5,2*nlayer+0.5],[0,length(dd)+1], ...
        'Color','w','linestyle','--','Linewidth',2);
    set(gca,'linewidth',2,'fontsize',20,'fontname','Times New Roman', ...
        'TickDir','out','TickLength',[0.007 0.01]);
    set(gcf,'unit','centimeters','position',[15,5,17,15]);

    text(-0.6,-1,'(c)','Color','black','fontsize',25,'fontname','Times New Roman');
    set(gcf,'PaperSize',[17 15]);
    saveas(gcf,'c.eps','psc2');
   
    %% 
    Jac = real(Jac);
    MYDATA_EXP = [MYDATA_EXP1;MYDATA_EXP2];
    Wd = diag(1./log(MYDATA_EXP(:,3))); 
    kn = 5;
    g = Jac'*(Wd'*Wd)*dd;
    for k = length(damp):-1:kn
        dm(:,k) = (Jac'*(Wd'*Wd)*Jac + damp(k)*(delta'*delta))\g;
        m(:,k) = mpre + dm(:,k);
        vs(:,k) = exp(real(m(1:nlayer,k)));
        vp(:,k) = exp(real(m(nlayer+1:2*nlayer,k)));
        v_rat(:,k) = vp(:,k)./vs(:,k);
        bos(:,k) = (0.5*v_rat(:,k).^2-1)./(v_rat(:,k).^2-1);

        [~,rd(k)] = guided_forward_vsvp(vs(:,k), vp(:,k), bos(:,k), rho, init_thk);
        fprintf('%%%%%%%%%% damp=%d %%%%%%%%%%', k);
        fprintf('\n');
    end

    %%
    rd(rd==0) = 100;
    [minrd, ind] = min(rd);
    
    mcur = m(:,ind);
    vscur = exp(real(mcur(1:nlayer)));
    vpcur = exp(real(mcur(nlayer+1:2*nlayer)));
    rdcur = minrd;
    varargout{1}(iter) = minrd;
    
    %%
    mresult = mcur; 
    vsresult = exp(real(mresult(1:nlayer)));
    vpresult = exp(real(mresult(nlayer+1:2*nlayer)));
    v_ratres = vpresult./vsresult;
    bosresult = (0.5*v_ratres.^2-1)./(v_ratres.^2-1);
    [~,varargout{1}(iter)] = guided_forward_vsvp(vsresult, vpresult, bosresult, rho, init_thk);
    
    fprintf('Interation:%d\n', iter);
    fprintf('|%.2f|', exp(mcur));
    fprintf('\n');
end
toc

%%
figure
for k = 1:MAXROOT
   H1=plot_vf(EXPDATA1{k}.freq, EXPDATA1{k}.vr, 200, 1800, '-', '[0.4 0.4 0.4]','.',3,25);
   H4=plot_vf(EXPDATA2{k}.freq, EXPDATA2{k}.vr, 200, 1800, '-', '[0.4 0.4 0.4]','.',3,25);
end
for k = 1:MAXROOT
   H2=plot_vf(INIDATA1{k}.freq, INIDATA1{k}.vr, 200, 1800, '-', '[0.2 0.7 0.5]','.',2,20); 
   H5=plot_vf(INIDATA2{k}.freq, INIDATA2{k}.vr, 200, 1800, '-', '[0.2 0.7 0.5]','.',2,20);
end
for k = 1:MAXROOT
   H3=plot_vf(THEDATA1{k}.freq, THEDATA1{k}.vr, 200, 1800, '-.', '[0.9 0.1 0.1]','^',2,3);
   H6=plot_vf(THEDATA2{k}.freq, THEDATA2{k}.vr, 200, 1800, '-.', '[0.9 0.1 0.1]','^',2,3);
end
legend([H1(1),H2(1),H3(1)],'True model','Initial model','Inversion result');
xlim([0,200]);
ylim([0,2000]);
set(gca,'YTick',200:400:1800);
set(gca,'XTick',0:50:200);
xlabel('Frequency (Hz)'); ylabel('Phase Velocity (m/s)');
set(gca,'linewidth',2,'fontsize',20);
grid on;
box on;
set(gca,'linewidth',2,'fontsize',20,'fontname','Times New Roman','TickDir','out','TickLength',[0.007 0.01]);
set(gcf,'unit','centimeters','position',[15,5,25,15]);

text(-27,2020,'(a)','Color','black','fontsize',25,'fontname','Times New Roman');
set(gcf,'PaperSize',[25 15]);
saveas(gcf,'a.eps','psc2');

%%
figure
hold on;
plotvsthk(true_vs,true_thk, 40, 0, 2500, '-','[0.4 0.4 0.4]','.',3,25);
plotvsthk(true_vp,true_thk, 40, 0, 2500, '-.','[0.4 0.4 0.4]','.',3,25);

plotvsthk(init_vs,init_thk, 40, 0, 2500, '-','[0.2 0.7 0.5]','.',2,20);
plotvsthk(init_vp,init_thk, 40, 0, 2500, '-.','[0.2 0.7 0.5]','.',2,20);

plotvsthk(vsresult,init_thk, 40, 0, 2500, '-', '[0.9 0.1 0.1]','^',2,3);
plotvsthk(vpresult,init_thk, 40, 0, 2500, '-.', '[0.9 0.1 0.1]','^',2,3);

xlabel('Depth (m)'); ylabel('Velocity (m/s)');
set(gca,'linewidth',2,'fontsize',20);
legend('True Vs','True Vp','Initial Vs','Initial Vp','Inverted Vs','Inverted Vp');
view(90,90);
grid on;
box on;
set(gca,'linewidth',2,'fontsize',20,'fontname','Times New Roman','TickDir','out','TickLength',[0.007 0.01]);
set(gcf,'unit','centimeters','position',[15,1,15,25]);

text(-0.9,-400,'(d)','Color','black','fontsize',25,'fontname','Times New Roman');
set(gcf,'PaperSize',[15 25]);
saveas(gcf,'d.eps','psc2');




