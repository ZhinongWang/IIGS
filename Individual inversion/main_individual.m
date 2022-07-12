clc; clear; close all;
global MYDATA EXPDATA THEDATA MAXROOT MYDATA_EXP MYDATA_THE Wd INIDATA init_vp

%% Read data
R0 = load("R0.mat");
R0 = R0.R0;
R1 = load("R1.mat");
R1 = R1.R1;
MYDATA = [R0;R1];
MYDATA_EXP = MYDATA;

%% Pre-Processing
preprocess;

%% Initial model
true_vs = [270 330 450 540];
true_vp = [900 1100 1300 1800];
true_thk = [5 10 5];

init_vs = [150 250 350 600];
init_vp = [700 900 1500 1900];
init_thk = [5 10 5];
rho = [1.6 1.6 1.8 2.0];
init_v_rat = init_vp./init_vs;
init_bos = (0.5*init_v_rat.^2-1)./(init_v_rat.^2-1);

%% Initial parameter
dv = 1.0;      
rdc = 0.01;    
maxrun = 20;  
init_m = [init_vs];

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
vpcur = init_vp';

for iter = 1:10
    v_ratio = vpcur/vscur;
    boscur = (0.5.*v_ratio.^2-1)./(v_ratio.^2-1);
    [~,rdcur] = guided_forward_vsvp(vscur, init_vp, boscur, rho, init_thk);
    if iter==1
        for k = 1:MAXROOT
            INIDATA{k}.freq = THEDATA{k}.freq;
            INIDATA{k}.vr = THEDATA{k}.vr;
        end
    end
    
    mpre = mcur;
    rdpre = rdcur;
    
    %% 
    dd = log(MYDATA_EXP(:,3)) - log(MYDATA_THE(:,3));
    dd(isinf(dd)) = 0.;
    Jac = guided_jacobian_vs(mcur, boscur, rho, init_thk);
    
    %%
    Jac_abs = abs(Jac);
    Jac_real = real(Jac);
    Jac_imag = imag(Jac);

    figure
    imagesc(Jac_real);   
    colormap(turbo);
    caxis([0 1.2]);
    colorbar;
    for j=1:MAXROOT
        line([0,length(mcur)+1],[EXPDATA{j}.ind(end)+0.5, EXPDATA{j}.ind(end)+0.5], ...
            'Color','w','linestyle','--','Linewidth',2);
    end
    line([nlayer+0.5,nlayer+0.5],[0,length(dd)+1], ...
        'Color','w','linestyle','--','Linewidth',2);
    line([2*nlayer+0.5,2*nlayer+0.5],[0,length(dd)+1], ...
        'Color','w','linestyle','--','Linewidth',2);
    set(gca,'linewidth',2,'fontsize',20,'fontname','Times New Roman','TickDir','out','TickLength',[0.007 0.01]);
    set(gcf,'unit','centimeters','position',[15,5,17,15]);
    
    text(-0.1,0,'(c)','Color','black','fontsize',25,'fontname','Times New Roman');
    set(gcf,'PaperSize',[17 15]);
    saveas(gcf,'c.eps','psc2');
   
    %%
    Jac = Jac_real;
    kn = 5;
    g = Jac'*(Wd'*Wd)*dd;
    for k = length(damp):-1:kn
        dm(:,k) = (Jac'*(Wd'*Wd)*Jac + damp(k)*(delta'*delta))\g;
        m(:,k) = mpre + dm(:,k);
        vs(:,k) = exp(real(m(1:nlayer,k)));
        vp(:,k) = init_vp;
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
    vpcur = init_vp';
    rdcur = minrd;
    varargout{1}(iter) = minrd;
    
    %% 
    mresult = mcur; 
    vsresult = exp(real(mresult(1:nlayer)));
    vpresult = init_vp;
    v_ratres = vpresult./vsresult;
    bosresult = (0.5*v_ratres.^2-1)./(v_ratres.^2-1);
    [~,varargout{1}(iter)] = guided_forward_vsvp(vsresult, vpresult, bosresult, rho, init_thk);
    
    fprintf('Iteration:%d\n', iter);
    fprintf('|%.2f|', exp(mcur));
    fprintf('\n');
end
toc
%%
figure
for k = 1:MAXROOT
   H1=plot_vf(EXPDATA{k}.freq, EXPDATA{k}.vr, 200, 1800, '-', '[0.4 0.4 0.4]','.',3,25); 
   H2=plot_vf(INIDATA{k}.freq, INIDATA{k}.vr, 200, 1800, '-', '[0.2 0.7 0.5]','.',2,20);  
   H3=plot_vf(THEDATA{k}.freq, THEDATA{k}.vr, 200, 1800, '-.', '[0.1 0.1 0.5]','^',2,3);
end
legend([H1(1),H2(1),H3(1)],'Observed data','Initial model','Inversion result');
xlim([0,100]);
ylim([0,600]);
set(gca,'YTick',0:150:600);
set(gca,'XTick',0:20:100);
xlabel('Frequency (Hz)'); ylabel('Phase Velocity (m/s)');
set(gca,'linewidth',2,'fontsize',20);
grid on;
box on;
set(gca,'linewidth',2,'fontsize',20,'fontname','Times New Roman','TickDir','out','TickLength',[0.007 0.01]);
set(gcf,'unit','centimeters','position',[15,5,17,15]);

text(-18,620,'(a)','Color','black','fontsize',25,'fontname','Times New Roman');
set(gcf,'PaperSize',[17 15]);
saveas(gcf,'a.eps','psc2');

%%
figure
hold on;
plotvsthk(true_vs,true_thk, 40, 0, 2500, '-','[0.4 0.4 0.4]','.',3,25);
plotvsthk(true_vp,true_thk, 40, 0, 2500, '-.','[0.4 0.4 0.4]','.',3,25);

plotvsthk(init_vs,init_thk, 40, 0, 2500, '-','[0.2 0.7 0.5]','.',2,20);
plotvsthk(init_vp,init_thk, 40, 0, 2500, '-.','[0.2 0.7 0.5]','.',2,20);

plotvsthk(vsresult,init_thk, 40, 0, 2500, '-', '[0.1 0.1 0.5]','^',2,3);
plotvsthk(vpresult,init_thk, 40, 0, 2500, '-.', '[0.1 0.1 0.5]','^',2,3);

xlabel('Depth (m)'); ylabel('Velocity (m/s)');
set(gca,'linewidth',2,'fontsize',20);
legend('True Vs','True Vp','Initial Vs','Initial Vp','Inverted Vs','Inverted Vp');
view(90,90);
grid on;
box on;
set(gca,'linewidth',2,'fontsize',20,'fontname','Times New Roman','TickDir','out','TickLength',[0.007 0.01]);
set(gcf,'unit','centimeters','position',[15,1,15,25]);

text(-0.9,-400,'(e)','Color','black','fontsize',25,'fontname','Times New Roman');
set(gcf,'PaperSize',[15 25]);
saveas(gcf,'e.eps','psc2');



