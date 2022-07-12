clc;
clear;
close all;
%%
nx = 4;
dx = 50;
dt = 0.001;

sis = load('Sis_4.mat');
sis4 = sis.sis4;

figure
imagesc(sis4);
caxis([-1 1]);
colormap(gray);
box on;
X = 1:100:860;
XX = X-1;
xticks(X);
xticklabels(XX);
Y = 1:1000:7001;
YY = (Y-1)*dt;
yticks(Y);
yticklabels(YY);
xlabel('Trace'); 
ylabel('Time (s)');

set(gca,'linewidth',2,'fontsize',20,'fontname','Times New Roman', ...
    'TickDir','out','TickLength',[0.007 0.01]);
set(gcf,'unit','centimeters','position',[1,1,40,15]);

set(gcf, 'PaperSize', [40.5 15.5]);
saveas(gcf,'1.eps','psc2');








