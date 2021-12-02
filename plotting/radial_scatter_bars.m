function radial_scatter_bars(f0,style,filled,color)
% Usage: radia_scatter_bars(f0,fig_handle=1)
%        plot scatter and error bars of f0 in figure fig_handle
% Columns of f0: [time mean min max]

if isempty(f0), return; end
% f0(:,2) = f0(:,2)/2;

if ~exist('style','var') || isempty(style), style = '*'; end
if ~exist('filled','var') || ~filled, filled = false; end
if ~exist('color','var') || isempty(color), color = 'k'; end

hold on
if filled
    scatter(f0(:,1),f0(:,2),color,style,'filled');
else
    plot(f0(:,1),f0(:,2),style,'color',color)
end
% plot(f0(:,1),f0(:,2),style,'color',color)

if size(f0,2)>2
    tickWidth = 0.05;
    for i=1:size(f0,1)
        barMean = f0(i,2);
        barAbove = f0(i,4) - f0(i,2);
        barBelow = f0(i,2) - f0(i,3);
        ebh = terrorbar(f0(i,1), barMean, barBelow, barAbove, tickWidth);
        set(ebh, 'color', color, 'linewidth', 1);
    end
end
% xlim([f0(1,1) f0(end,1)+1])