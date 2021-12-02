function radial_plot_single_run(ylims,DIRBASE)
global OUTBASE legends idx ts makeMovies R record_every
close all

% DIRBASE = 'H:\Dropbox\Dropbox\Research\radial_param_study\feedback on lambda\';

% time points to plot. Comment to plot at times from data
ts = [15:5:30 40:20:80 100];
% ts = [5:5:30 50 100];
% ts = 12:2:28;
% ts = [1:0.25:2 3:5];

plot_R_only_time_pts = 0; % only plot R(t) up to the latest data time point

% if exist('ylims','var')
%     p_ylim = ylims(1,:);
%     hoop_ylim = ylims(2,:);
%     radial_ylim = ylims(3,:);
%     sigma_rr_ylim = ylims(4,:);
%     sigma_tt_ylim = ylims(5,:);
% end
R_ylim = [0 600];
p_ylim = [];
hoop_ylim = [];
radial_ylim = [];
sigma_rr_ylim = [];
sigma_tt_ylim = [];
prv_ylim = [];
rel_hoop_ylim = [];
rel_radial_ylim = [];

makeMovies = 0;

OUTBASE = './';
if ~exist('DIRBASE','var'), DIRBASE='.'; end
load([DIRBASE '/solution.mat']);
if exist('solution_rhoc.mat','file')
    load solution_rhoc RHOC LAC
end
if ~exist([DIRBASE '/parameters.mat'], 'file')
    warning('can''t load parameters.mat. assuming records at integer frames');
    idx = ts;
else
    load([DIRBASE '/parameters.mat']);
    if ~exist('ts','var'), ts = f0(:,1); end
    if ~exist('record_every','var'), record_every = 1; end
    dataSize = (tspan(2)-tspan(1))/dt;
    dataFrameFreq = (tspan(2)-tspan(1))/(dataSize-1);
    tsFrameReq = min(ts(2:end)-ts(1:end-1));
    if tsFrameReq<dataFrameFreq-1e-5
        error(['insufficient recorded frequency. suggested numFrames = ' num2str((tspan(2)-tspan(1))/tsFrameReq)]);
    end
    idx = int32((ts-tspan(1))/dataFrameFreq+1)/1;
    idx = min(max(idx,2),numFrames-2);   % do not plot 1st and last frame
    if sum(idx>dataSize)>0, error('requested frame exceeds data size'); end
end

% R(t)
h = figure('visible','off');
plot(linspace(tspan(1),tspan(2),length(R)),R,'linewidth',2);
xl = xlim;
% f0(:,2) = f0(:,2)/2;
radial_scatter_bars(f0); xlim(xl);
% if exist('R_ylim','var') && ~isempty(R_ylim), ylim(R_ylim); end
title('Tumor radius (\mum)'); xlabel('Time');
if plot_R_only_time_pts, xlim([xl(1) f0(end,1)]); end
print(h,'-dpng','-r300',[OUTBASE 'R vs T']); close(h);
% return

% make legends
if ~makeMovies
    legends = cell(1,length(ts));
    for i=1:length(ts), legends{i}=['T=' num2str(ts(i))]; end
    h = figure('visible','off'); hold on
    for i=1:length(ts), plot(1,1,'linewidth',max(5-0.75*i, 0.75)); end
    axis off
    legend(legends,'interpreter','latex','location','northwest');
    print(h,'-dpng','-r300',[OUTBASE 'legends']); close(h);
end

% r y p v R radial hoop VT RHO int_rhoc YR C B LA LAMBDA TMP
if ~exist('del_radial','var')
    do_plot(r,C,'c(r,t)');
    do_plot(r,y,'y(r,t)');
    do_plot(r,v,'v(r,t)',prv_ylim,'plot(r,zeros(size(r)),''k'',''linewidth'',0.5)');
    do_plot(r,p,'p(r,t)',p_ylim);
% %     
    do_plot(r,radial,'\sigma_{rr} - p',radial_ylim,'plot(r,zeros(size(r)),''k'',''linewidth'',0.5)');
    do_plot(r,hoop,'\sigma_{\theta\theta} - p',hoop_ylim,'plot(r,zeros(size(r)),''k'',''linewidth'',0.5)');
    do_plot(r,radial+p,'\sigma_{rr}',sigma_rr_ylim);
    do_plot(r,hoop+p,'\sigma_{\theta\theta}',sigma_tt_ylim);
% 
%     do_plot(r,-PD,'-p^D');
%     do_plot(r,(hoop+PD),'hoop+p^D',rel_hoop_ylim);
%     do_plot(r,(radial+PD),'radial+p^D',rel_radial_ylim);
    
%     do_plot(r,PRA./R(1:length(R)-1)','p_r over R');
else
%     do_plot(r,del_radial,'\Delta(\sigma_{rr} - p)',radial_ylim);
%     do_plot(r,del_hoop,'\Delta(\sigma_{\theta\theta} - p)',hoop_ylim);
end
% 
% do_plot(r,lamBs,'\lambda_B');
do_plot(r,LAMBDA,'\lambda');
do_plot(r,LAMBDA.*C,'\lambda c');
do_plot(r,LA,'\lambda_A');
do_plot(r,LAMBDA.*C-LA,'\lambda c-\lambda_A',[],'plot(r,zeros(size(r)),''k'',''linewidth'',0.5)');

radial_plot_rhoc(ts);

% do_plot(r,TMP,'st');
% do_plot(r,VT,'rescaled v')
% do_plot(r,YR,'y_r');
% do_plot(r,(YR./y.*r).^(2/3),'x');

% lambdaA = lambdaA_base + fcA * (lambdaA_A*abs(hoop-s0cA)).^nlamA.*(hoop<s0cA);
% do_plot(r,lambdaA,'\lambda_A');
% do_plot(r,Lbase./(1+(L_A*(hoop-s0c)).^nL.*(hoop<s0c)),'L');
% do_plot(r,lambda_base./(1+(lambda_A*(hoop-s0c)).^nlam.*(hoop<s0c)),'\lambda');
% do_plot(r,C.*lambda_base./(1+(lambda_A*(hoop-s0c)).^nlam.*(hoop<s0c)),'\lambda c');
end

function do_plot(r,data,name,ylim_,postPlotCmd)
global OUTBASE legends idx makeMovies R record_every

if makeMovies, make_movie(r,data,name); return; end

h = figure('visible','off'); hold on
for i=1:length(idx)
    this = max(idx(i),2);
    plot(r*R(this),data(:,this),'linewidth',max(5-0.75*i, 0.75));
end
if exist('ylim_','var') && ~isempty(ylim_), ylim(ylim_); end
ylabel(name);
xlabel('R');
if exist('postPlotCmd','var'), eval(postPlotCmd); end
% legend(legends,'location','northeast');
name = strrep(strrep(name,'\',''),'/',' over ');
print(h,'-dpng',[OUTBASE name '.png'],'-r300'); close(h);
end

function make_movie(r,data,name)
global OUTBASE idx ts
if isempty(ts), return; end

disp(['Making movies, ' name]);
if ispc, png_dir = 'I:/radial_png';
else, png_dir = '~/w1/radial/radial_png'; end

warning off MATLAB:MKDIR:DirectoryExists
mkdir(png_dir)
warning on MATLAB:MKDIR:DirectoryExists
delete([png_dir '/*.png']);

% calculate ylim automatically. 3-7-17
data_ = data(:,2:end-1);
ylim_ = [min(min(data_)) max(max(data_))];

for i=1:length(idx)
    h = figure('visible','off');
    plot(r,data(:,idx(i)),'linewidth',2);
    ylim(ylim_);
    title([name ', T=' num2str(ts(i))]);
%     xlabel('r''');
    print(h,'-dpng',[png_dir '/' num2str(ts(i)) '.png'],'-r250');
    close(h);
end

mov = [OUTBASE '_' strrep(name,'\','') '.avi'];
writerObj = VideoWriter(mov);
writerObj.FrameRate = 2;
open(writerObj);
for i=1:length(idx)
    filename = [png_dir '/' num2str(ts(i)) '.png'];
    if ~exist(filename,'file'), continue; end
    thisimage = imread(filename);
    writeVideo(writerObj, thisimage);
end
close(writerObj);
ffmpeg(mov);
end

% ts = [8:2:20 30 40]; % Helminger 1.0%
% ts = 8:8:80;
% ts = 4:4:20;
% ts = [4:4:24 30 40 60 80];
% ts = [20 29 30 31 35 40 50 80 100];
% ts = [5 10 20 30 40 60 100 200];
% ts = [15:5:30 40:20:80 100 ];
% ts = [1 5 10 20 30 40 200];
% ts = [8:4:24 30 40 60 80 100 150 200];
% ts = [100 150 200];
% ts = [4:0.2:5 6 8 10];
% ts = [30 50 80];
% ts = [1:5 7 9 12 15 20 30 40 100];

% ts = 1:0.25:3;
% ts = 1:10;
% ts = [1:0.25:2 3 4 5];
% ts = 4:10;
% ts = [3 5 7 9 15 20 40 60 80 100];