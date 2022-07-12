function [H]=plot_vf(freq, vr, vfirst, vlast, mstyle, color, mkr, lir, sir)
%   画频散曲线
%   vr -- 相速度或群速度频散曲线(每列对应一个模式)
%   freq -- 频率范围
%   style  color -- 指定曲线的风格和颜色

mn = size(vr,2); % 每列存放一个模式
for i = 1:mn
    ind = find(vr(:,i) > 0);
    H=plot(freq(ind), vr(ind,i),'Marker',mkr, 'LineStyle', mstyle, 'Color', color,'linewidth',lir,'MarkerSize',sir);
    hold on;
end
ylim([vfirst-0.1,vlast]); % 速度显示范围
xlabel('Frequency(Hz)');
ylabel('Phase Velocity(m/s)');
%hold off;
end

