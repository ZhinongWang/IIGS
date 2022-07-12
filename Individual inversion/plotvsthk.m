function fh = plotvsthk(vs,thk,depthest,v0,v1,style,clr,mkr,lir,sir)
% 画速度层厚度阶梯曲线
% fh--返回图形句柄 depthest--图形显示最大深度
% v0,v1--图形显示的速度范围
gcf;
vs = vs(:); thk = thk(:);
depth = cumsum(thk);
if nargin == 7
    fh = stairs([0;depth;depthest],[vs;vs(end)],'LineStyle',style,'Color',clr,'LineWidth',lir);
else
    fh = stairs([0;depth;depthest],[vs;vs(end)],'LineStyle',style,'Color',clr,'Marker',mkr,'LineWidth',lir,'MarkerSize',sir);
end

xlim([0 depthest]);
ylim([v0 v1]);
end