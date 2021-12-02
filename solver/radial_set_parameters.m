function radial_set_parameters(param)

lambda_base = param(1);
lambdaA_base = param(2);
Lbase = param(3);
lambda_mr = param(4);
% c_lam = param(4); % see get_c
% L_A = param(4);
lambda_A = param(5);
pBar = param(6);
s0c = param(7);
nL = param(8);
beta_base = param(9);
% beta_w = param(10);
c_lamB = param(10); % see get_c
lambdaA_A = param(11); % delta_lamA
fcA = param(12); % A (sensitivity)
cH = param(13);
tumorID = param(14);
mu = param(15); % mooney-rivlin
cH2 = param(16); % mooney-rivlin, host
gamma_B = param(17); % feedback on lam_B. see get_c

gLamMns = param(18); %negative feedback strength on lambda
nLamMns = param(19); %exponent associated with negative feedback for lambda
lambda_max = param(20); %max effect of positive feedback on lambda
gLamPls = param(21);%positive feedback strength on lambda
nLamPls = param(22);%exponent associated with posotive feedback for lambda

cB = 1; % c at boundary and in vessels

cT = 1; gamma_ = 0;
% beta_base = 1; beta_w = 0; % w=0: mitosis only; w=1: apoptosis only
scale_v = 0.001;
with_G_incompatibility = 1;

L_A = 0;
beta_w = 0;

% Nr = 1e4; dt = 5e-5; record_every = 200; % record every this many time steps
Nr = 5e2; dt = 1e-2; record_every = 1; % record every this many time steps
% Nr = 200; dt = 1d-2;

% if beta_base<0.1, Nr = 5e2; dt = 0.02; end
% if lambda_base>1.5, dt = 0.01; end

switch tumorID
    case -2 % parameter study of Morgan. load Free data regardless of pBar
        load morgan_fig1c_free1s.mat f0
        f0 = f0(2:end,:);
        tspan = [f0(1,1) 100];
    case -1
        isFree = pBar==0;
        if isFree
            load morgan_fig1c_free1s.mat f0
            f0 = f0(2:end,:);
        else
            load morgan_fig1c_pressures.mat f0
            f0 = f0(2:end,:);
        end
%         tspan = [0 1];
        tspan = [f0(1,1) 100];
    case 0
        load helminger_fig1a.mat f0
%         tspan = [f0(1,1) 200];
        tspan = [f0(1,1) 200];
    case {1,7}
        load helminger_fig1a.mat f07
        f0 = f07;%(2:end,:);
        tspan = [f0(1,1) 200];
    case 2
        load helminger_fig1a.mat f10
        f0 = f10;
        tspan = [f0(1,1) 200];
    case 3
        load helminger_fig1a.mat f03
        f0 = f03(1:end,:);
        tspan = [f0(1,1) 200];
    case 5
        load helminger_fig1a.mat f05
        f0 = f05;%(2:end,:);
        tspan = [f0(1,1) 200];
    case 8
        load helminger_fig1a.mat f08
        f0 = f08;%(2:end,:);
        tspan = [f0(1,1) 200];
    case 9
        load helminger_fig1a.mat f09
        f0 = f09;%(2:end,:);
        tspan = [f0(1,1) 200];
    case 10
        load helminger_fig1b.mat f0
        f0 = f0(2:end,:);
        tspan = [f0(1,1) 150];
    case 11
        load helminger_fig1b.mat f07
        f0 = f07(3:end,:);
        tspan = [f0(1,1) 150];
    case 12
        load helminger_fig1b.mat f10
        f0 = f10(:,:);
        tspan = [f0(1,1) 150];
    case 20
        load prl_fig1as.mat f0
        tspan = [f0(1,1) 100];
    case 21 % Morgan pressure release
        load prl_fig1a_500pa_releases.mat f0
        f0 = f0(2:end,:); % skip first data point
        tspan = [f0(1,1) 100];
    case 221 % Morgan pressure release, pre-release data
        load prl_fig1a_500pa_releases.mat f0
        f0 = f0(2:end-1,:);
        tspan = [f0(1,1) 100];
end

if exist('radial_grid_search.lock','file')
    tspan = [f0(1,1) f0(end,1)];
end

R0 = f0(1,2);
% R0 = 50;

% tspan = [0 20]; dt = 2e-2; R0 = 175; % PRL Fig 1a free
% tspan = [3 20]; dt = 1e-2; R0 = 173.5293; % PRL Fig 1a 500pa
% tspan = [0 300];  dt = 1e-2; R0 = 75;

% Lbase = 100; lambda_base = 1; lambdaA_base = 0.5;
% Lbase = 210; lambda_base = 0.8; lambdaA_base = 0.5; % PRL Fig 1a
% Lbase = 170; lambda_base = 0.7; lambdaA_base = 0.5; % morgan free

% cH = 0; %pBar = 0;
% lambdaA_A = 0; s0cA = 0;
% lambda_A = 10; L_A = 1e4;
% s0c = -.35;
s0cL = s0c;
s0cA = s0c;
% nL = 2;
nlam = nL; nlamA = nL;
% fcA = 0; % fold change of feedback on lambda_A

lambdaC = lambda_base;

match_R_only = 0;

disp_progress = ~exist('radial_grid_search.lock','file');
numFiguresSamePlot = 8;

% numFrames = max(numFiguresSamePlot, tspan(2)-tspan(1)); % only store these many time points beyond tspan(1)
% numFrames = max(numFiguresSamePlot, 800);
% if numFrames > (tspan(2)-tspan(1))/dt, error('Requested recording frequency too high. Reduce numFrame or dt.'); end
numFrames = (tspan(2)-tspan(1))/dt;

if ~exist('f0','var'), f0 = []; end

if ~exist('record_every','var'), record_every = 1; end
frame2id = @(x) int32((x-tspan(1))/dt/record_every);

save('parameters.mat', ...
    'param','disp_progress','match_R_only','numFrames','numFiguresSamePlot','f0',...
    'R0','cT','cH','gamma_','dt','tspan','Nr','beta_base','beta_w','Lbase','with_G_incompatibility',...
    'lambda_base','lambdaA_base','lambdaC','scale_v','pBar',...
    'lambdaA_A','s0cA','nlamA','lambda_A','s0c','nlam','L_A','s0cL','nL','fcA',...
    'mu','cH2','tumorID','lambda_mr','c_lamB','gamma_B','cB',...
    'gLamMns','nLamMns','lambda_max','gLamPls','nLamPls','record_every','frame2id');

end
