function solve_rhoc_fast(param,do_not_save)
% read existing solutions and solve rhoc equation only

load parameters fcA nlamA lambdaA_base lambdaA_A dt tspan record_every ...
    Nr numFrames lambda_A lambda_base nlam lambda_mr ...
    gLamMns nLamMns lambda_max gLamPls nLamPls
    
load solution v r R hoop C 

if exist('param','var') && ~isempty(param)
    Drho = param(1); gamma_ac = param(2); nrho = param(3);
else
    Drho = 0.;
    gamma_ac = 0.1; 
    nrho = 2;
end

% restart = 200;

Nt = (tspan(2)-tspan(1))/dt;
dx = r(2)-r(1);

if ~exist('do_not_save','var'), disp_progress = 10; else, disp_progress = 0; end
if exist('radial_grid_search.lock','file'), disp_progress = 0; rhoc_warned = true; end
if ~exist('record_every','var'), record_every = 1; end

disp_progress = 0;

LC = zeros(Nr, int32(numFrames/record_every+1));
LAC = zeros(Nr, int32(numFrames/record_every+1));
RHOC = zeros(Nr, int32(numFrames/record_every+1));
vdr = zeros(Nr, int32(numFrames/record_every+1));
VT = zeros(Nr, int32(numFrames/record_every+1));

rhoc = ones(Nr,1);
RHOC(:,1) = rhoc;
    
for n=1:Nt
    currTime = tspan(1)+dt*(n-1);
    
    if disp_progress && mod(currTime-tspan(1),disp_progress)==0, disp(['T=' num2str(currTime)]); end
    
    vt = (v(:,n)-r*v(end,n))/R(n);
    vt = vt*1; % change the effect of v*div(rho_c)
    
    fs = (abs(hoop(:,n)).*(hoop(:,n)<0)).^nlamA;
    lambdaA = lambdaA_base + lambdaA_A * fcA*fs./(1+fcA*fs);
%     fs = (lambda_A*abs(hoop(:,n))).^nlam.*(hoop(:,n)<0);
%     lambda = lambda_base ./ (1 + fs);
%     fs = (abs(hoop(:,n)).*(hoop(:,n)>0)).^nlam;
%     lambda = lambda_base;% + lambda_mr * lambda_A*fs./(1+lambda_A*fs);
        
    lambda = lambda_base./(1+gLamMns*(abs(hoop(:,n)).*(hoop(:,n)<0)).^nLamMns)...
        +lambda_max*(hoop(:,n) >0).*(gLamPls*(abs(hoop(:,n))).^nLamPls)./...
        (1+gLamPls*(abs(hoop(:,n))).^nLamPls);

    lambdac = lambda_base ./ (1 + gamma_ac*abs(rhoc-1).^nrho);
    lambdaAc = lambdaA_base;
%     thres_c = -0.66;
%     fsc = (abs(hoop(:,n)).*(hoop(:,n)<thres_c)).^nlamA;
%     lambdaAc = lambdaA_base + lambdaA_A * fcA*fsc./(1+fcA*fsc);
    
    if exist('restart','var') && currTime<restart, rhoc = ones(Nr,1); continue; end
    
    tmp1 = Drho/R(n)^2/(2*dx*dx);
    rhoc_new = ones(Nr,1);
    
    rhs = (lambdac-lambda).*C(:,n) + (lambdaA-lambdaAc);
    vp = max(vt, 0); vn = min(vt, 0);
    
    % at boundary, vt=0, homogeneous Neumann
    upw = - tmp1*(rhoc(Nr-1)*(rhoc(Nr-1)-rhoc(Nr))-rhoc(Nr-1)*(rhoc(Nr-1)-rhoc(Nr-2)));
    rhoc_new(end) = (rhoc(end)-dt*upw)/(1-dt*rhs(end));

    % homogeneous Neumann at origin, rhoc(2)-rhoc(0)=0
    upw = (vp(1)*(rhoc(1)-rhoc(2)) + vn(1)*(rhoc(2)-rhoc(1)))/dx ...
        - tmp1*(rhoc(2)*(rhoc(3)-rhoc(2))-rhoc(2)*(rhoc(1)-rhoc(2)));
    rhoc_new(1) = (rhoc(1)-dt*upw) / (1-dt*rhs(1));
    
    for i=2:Nr-2
        upw = (vp(i)*(rhoc(i)-rhoc(i-1)) + vn(i)*(rhoc(i+1)-rhoc(i)))/dx ...
            - tmp1*(rhoc(i+1)*(rhoc(i+2)-rhoc(i+1))-rhoc(i-1)*(rhoc(i)-rhoc(i-1))) ...
            - 2/R(n)^2/r(i)*Drho*rhoc(i)*(rhoc(i+1)-rhoc(i-1))/(2*dx);
        rhoc_new(i) = (rhoc(i)-dt*upw) / (1-dt*rhs(i));
    end
    
    upw = (vp(Nr-1)*(rhoc(Nr-1)-rhoc(Nr-2)) + vn(Nr-1)*(rhoc(Nr)-rhoc(Nr-1)))/dx ...
        - tmp1*(rhoc(Nr)*(rhoc(Nr)-rhoc(Nr-1))-rhoc(Nr-2)*(rhoc(Nr-1)-rhoc(Nr-2))) ...
        - 2/R(n)^2/r(Nr-1)*Drho*rhoc(Nr-1)*(rhoc(Nr)-rhoc(Nr-2))/(2*dx);
    rhoc_new(Nr-1) = (rhoc(Nr-1)-dt*upw) / (1-dt*rhs(Nr-1));
    
    rhoc = rhoc_new;
    
    if ~exist('rhoc_warned','var') && max(abs(rhoc))>1e3
        warning(['rhoc not converging, T=' num2str(currTime)]); 
        rhoc_warned = true;
    end
    
    if mod(n,record_every)==0
        nn = n/record_every;
        RHOC(:,nn) = rhoc;
        LC(:,nn) = lambdac;
        LAC(:,nn) = lambdaAc;
        vdr(:,nn) = vt.*[0; rhoc(2:end)-rhoc(1:end-1)]/dx;
        VT(:,nn) = vt;
    end
end

if (~exist('do_not_save','var') || ~do_not_save) 
    if ~exist('radial_grid_search.lock','file')
        save solution_rhoc.mat RHOC LAC LC vdr VT
%         radial_plot_rhoc;
    else
        save solution_rhoc.mat RHOC 
    end
end