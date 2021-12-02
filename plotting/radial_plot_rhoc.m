function radial_plot_rhoc(ts)
global OUTBASE legends idx makeMovies R record_every
close all

% time points to plot
if nargin<1
%     ts = [8:4:24 30 40 60 100 150 200];
%     ts = [15:5:30 40:20:80 100 ];
%     
%     ts = [4:0.5:6 7 8];
%     ts = [9:0.5:11 12:15];
%     ts = [5:0.5:7 8 9];
%     
%     ts = [30 50 80];
%     ts = 10:10:80;
%     ts = [10:16];
%     ts = [200 201 202:2:210 220 250 300 400];

%     ts = [1:2:9 12 15 20 30 40 100];
%     ts = 0:0.1:1;
%     ts = [1:0.25:2 3 4 5];
%     ts = [3 5 7 9 15 20 40 60 80 100];
ts = [1:0.25:2 3:5];
end

if ~exist('solution_rhoc.mat','file'), warning('solution_rhoc.mat not found'); return; end

load solution_rhoc RHOC LAC LC vdr
load solution R r LA hoop LAMBDA C

makeMovies = 0;

OUTBASE = './';
if ~exist('DIRBASE','var'), DIRBASE='.'; end
if ~exist([DIRBASE '/parameters.mat'], 'file')
    warning('can''t load parameters.mat. assuming records at integer frames');
    idx = ts;
else
    load([DIRBASE '/parameters.mat']);
    if ~exist('ts','var'), ts = f0(:,1); end
    if ~exist('record_every','var') || isempty(record_every), record_every = 1; end
    dataSize = (tspan(2)-tspan(1))/dt;
    dataFrameFreq = (tspan(2)-tspan(1))/(dataSize-1);
    tsFrameReq = min(ts(2:end)-ts(1:end-1));
    if tsFrameReq<dataFrameFreq-1e-5
        error(['insufficient recorded frequency. suggested numFrames = ' num2str((tspan(2)-tspan(1))/tsFrameReq)]);
    end
    idx = int32((ts-tspan(1))/dataFrameFreq/record_every+1)/1;
%     idx = min(max(idx,2),numFrames-2);   % do not plot 1st and last frame
    if sum(idx>dataSize)>0, error('requested frame exceeds data size'); end
end

% % integral of rhoc
% intrhoc = zeros(1,size(RHOC,2)-1);
% for i=1:size(RHOC,2)-1
%     intrhoc(i) = simps(R(i*record_every)*r,(R(i*record_every)*r).^2.*RHOC(:,i));
% end
% save solution_rhoc intrhoc -append
% h = figure('visible','off');
% plot(linspace(tspan(1),tspan(2),length(intrhoc)),intrhoc);
% xlabel('T')
% print(h,'-dpng','-r300',[OUTBASE 'integral of rhoc']); close(h);

% % integral of RHS1
% intrhs1 = zeros(1,size(RHOC,2)-1);
% RHS1 = (LC.*C-LAC).*RHOC;
% for i=1:size(RHOC,2)-1
%     intrhs1(i) = simps(r,r.^2.*RHS1(:,i));
% end
% save solution_rhoc intrhs1 -append
% h = figure('visible','off');
% plot(linspace(tspan(1),tspan(2),length(intrhs1)),intrhs1);
% xlabel('T')
% print(h,'-dpng','-r300',[OUTBASE 'integral of RHS1']); close(h);

% make legends
legends = cell(1,length(ts));
for i=1:length(ts), legends{i}=['T=' num2str(ts(i))]; end
h = figure('visible','off'); hold on
for i=1:length(ts), plot(1,1,'linewidth',max(5-0.75*i, 0.75)); end
axis off
legend(legends,'interpreter','latex','location','northeast');
print(h,'-dpng','-r300',[OUTBASE 'legends']); close(h);

if size(LA,2)>size(LAC,2), LA=LA(:,1:size(LAC,2)); end

% do_plot(r,vdr,'v drho dr');
do_plot(r,(LC-LAMBDA).*C,'(\lambda_c-\lambda)n');
do_plot(r,LA-LAC,'\lambda_A-\lambda_{A,c}');
do_plot(r,(LAMBDA.*C-LA).*RHOC,'(\lambda n-\lambda_A)\rho');
do_plot(r,(LC.*C-LAC).*RHOC,'(\lambda_c n-\lambda_{A,c})\rho');
do_plot(r,((LC-LAMBDA).*C+(LA-LAC)).*RHOC-vdr,'RHS rho-advection',[],'plot(xlim,xlim*0,''k'',''linewidth'',0.5)');

% do_plot(r,C,'nutrient');
% do_plot(r,v,'v');
% do_plot(r,VT,'rescaled v');
do_plot(r,LA,'\lambda_A');
do_plot(r,LC,'\lambda_c');
do_plot(r,LAC,'\lambda_{A,c}');
do_plot(r,((LC-LAMBDA).*C+(LA-LAC)).*RHOC,'RHS rho',[],'plot(xlim,xlim*0,''k'',''linewidth'',0.5)');
do_plot(r,((LC-LAMBDA).*C+(LA-LAC)),'RHS',[],'plot(xlim,xlim*0,''k'',''linewidth'',0.5)');
% do_plot(r,((LC-LAMBDA)+(LA-LAC)),'RHS without n',[],'plot(xlim,xlim*0,''k'',''linewidth'',0.5)');
do_plot(r,RHOC,'\rho_c');
do_plot(r,hoop,'\sigma_{\theta\theta} - p');
end

function do_plot(r,data,name,ylim_,postPlotCmd)
global OUTBASE legends idx makeMovies R record_every

if makeMovies, make_movie(r,data,name); return; end

h = figure('visible','off'); hold on
for i=1:length(idx)
    this = max(idx(i),2);
    Ridx = this*record_every;
    plot(r*R(Ridx),data(:,this),'linewidth',max(5-0.75*i, 0.75));
end
if exist('ylim_','var') && ~isempty(ylim_), ylim(ylim_); end
ylabel(name);
xlabel('R');
if exist('postPlotCmd','var'), eval(postPlotCmd); end
% legend(legends,'location','northeast');
name = strrep(strrep(name,'\',''),'/',' over ');
print(h,'-dpng',[OUTBASE name '.png'],'-r300'); close(h);
end
