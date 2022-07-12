function guidedphase(FREQ,VS,VP,BOS,RHO,THK)
global EXPDATA THEDATA MAXROOT MYDATA_THE MYDATA_EXP VRN

VS = VS(:);
VP = VP(:);
BOS = BOS(:);
RHO = RHO(:);
THK = THK(:);
cn = length(VP);

%%
FN = length(FREQ);
% vrt = zeros(FN,MODE_NUMBER);

DVR = 0.1;                
vrmin = min(VS)*0.7;
vrmax = max(VS);
VR = vrmin:DVR:vrmax;
vrn = length(VR);

DVI = 1;
vimin = -10;
VI0 = -(vimin/DVI)+1;
vimax = 10;
VI = vimin:DVI:vimax;
vin = length(VI);

%%
cc = 1;
for fi = 1:FN
    fprintf('Freqency=%f\n',FREQ(fi));    
    parfor vri=1:vrn
        for vii=1:vin
            V = VR(vri)+1i*VI(vii);
            dpvalue = FastDelta1(VP,VS,RHO,THK,FREQ(fi),V);             
            dp(vri,vii) = dpvalue;
        end
    end
    dpa = abs(dp);
    dpal = log(dpa);
    %%
    isMin = islocalmin(dpa,1) & islocalmin(dpa,2);
    w = 2*pi*FREQ(fi);
    [row,col,v] = find(isMin);
    lr = length(row);
    for lri = 1:lr
        vv = VR(row(lri))+1i*VI(col(lri));
        kkk(cc) = w/vv;
        vv = w/real(kkk(cc));        
        if vv<max(VP)
            f1(cc) = FREQ(fi);
            v1(cc) = w/real(kkk(cc));
            a1(cc) = -imag(kkk(cc));
            x1(cc) = col(lri);
            y1(cc) = row(lri);            
            z1(cc) = dpa(row(lri),col(lri));
            cc = cc+1;
        end
    end
end

%%
cnum = cc-1;
% for i=1:cnum
%     if a1(i)==0
%         f1(i) = f1(i+1);
%         v1(i) = v1(i+1);
%         a1(i) = a1(i+1);
%         x1(i) = x1(i+1);
%         y1(i) = y1(i+1);
%         z1(i) = z1(i+1);
%     end
% end
toc;
%%
% figure
% % scatter(f1,v1,20,'r','filled');
% scatter(f1,v1,20,a1,'filled');
% % scatter(f1,v1,20,f1,'filled');
% xlim([0,200]);
% ylim([200,1800]);
% xlabel('Frequency (Hz)'); ylabel('Phase Velocity (m/s)');
% set(gca,'linewidth',2,'fontsize',20);
% grid on;
% box on;
% set(gca,'linewidth',2,'fontsize',20,'fontname','Arial','TickDir','out','TickLength',[0.007 0.01]);
% set(gcf,'unit','centimeters','position',[15,5,17,15]);

%%
% figure()
% imagesc(dpal);
% hold on;
% % line([0,vin],[241.5,241.5],'Color','w','linestyle','--','Linewidth',2);
% % line([0,vin],[1493,1493],'Color','w','linestyle','--','Linewidth',2);
% line([VI0,VI0],[0,vrn],'Color','w','linestyle','--','Linewidth',2);
% colormap jet;
% hold on;
% scatter(x1,y1,20,a1,'filled');
% % scatter(x,y,20,f1,'filled');
% axis tight;
% box on;
% view(0,-90);
% % caxis([0 100]);
% xlabel('Imag Phase Velocity (m/s)'); 
% ylabel('Real Phase Velocity (m/s)');
% title('dpal');
% set(gca,'linewidth',2,'fontsize',25,'fontname','Times New Roman');
% set(gcf,'unit','centimeters','position',[15,5,17,15]);

%% Data arrangement
for m = 1:MAXROOT
    lF = length(EXPDATA{m}.freq);
    n = 1;
    for i = 1:lF   %该模式观测数据频率
        k = 1; 
        for j = 1:cnum  %所有正演数据
            if EXPDATA{m}.freq(i)==f1(j)
                vv(k) = v1(j);
                ff(k) = f1(j);
                k = k+1;
            end
        end
        
        if exist('vv','var')
            vvs = unique(vv);
            vrt(n) = vvs(m);
            n = n+1;
        else
            vrt(n) = vrt(n-1);
            n = n+1;
        end

        clear vv vvs;
    end

    THEDATA{m}.vr = vrt(:);
    MYDATA_THE(EXPDATA{m}.ind(1):EXPDATA{m}.ind(end),3) = THEDATA{m}.vr;
    clear vrt;
end

%%
figure;
for k = 1:MAXROOT
   H1=plot_vf(EXPDATA{k}.freq, EXPDATA{k}.vr, 200, 1800, '-', 'b');
   H2=plot_vf(THEDATA{k}.freq, THEDATA{k}.vr, 200, 1800, '-', 'r');
end
legend([H1(1),H2(1)],'Observed data','Inversion result');
xlim([0,200]);
ylim([200,1800]);
set(gca,'YTick',200:400:1800);
set(gca,'XTick',0:50:200);
xlabel('Frequency (Hz)'); ylabel('Phase Velocity (m/s)');
set(gca,'linewidth',2,'fontsize',20);
grid on;
box on;
set(gca,'linewidth',2,'fontsize',20,'fontname','Arial','TickDir','out','TickLength',[0.007 0.01]);
set(gcf,'unit','centimeters','position',[15,5,17,15]);





















