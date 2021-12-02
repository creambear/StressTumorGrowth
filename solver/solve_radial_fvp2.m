function [r,Y,P,V,R,radial,hoop,VT,YR,C,B,LA,LAMBDA,TMP,lamBs,VA,PRA,CA,PD,residual] = solve_radial_fvp2(~)
if nargin==0, radial_time_evolution; return; end % call radial_time_evolution when directly executed

load('parameters.mat', ...
    'disp_progress','numFrames','numFiguresSamePlot','f0','record_every',...
    'R0','cT','cH','gamma_','dt','tspan','Nr','beta_base','Lbase','with_G_incompatibility',...
    'lambda_base','lambdaA_base','lambdaC','scale_v','pBar',...
    'lambdaA_A','s0cA','nlamA','lambda_A','s0c','nlam','s0cL','nL','fcA',...
    'lambda_mr','c_lamB','gamma_B','cB','gLamMns','nLamMns','lambda_max','gLamPls','nLamPls');%,'mu','cH2');

if ~exist('mu','var'), mu = 0; end
if ~exist('cH2','var'), cH2 = 0; end

R0_ref = 75;
% R0_ref = R0;

if ~exist('record_every','var'), record_every = 1; end

% Helminger Fig 1b pressure release time. 0 to disable.
if ~ismember('match_fig1b', who('-file', 'parameters.mat'))
    match_fig1b = 0; 
else
    load('parameters.mat','match_fig1b');
    if match_fig1b>0 && ~exist('radial_grid_search.lock','file')
%         disp(['Pressure release at T=' num2str(match_fig1b)]);
        R0_ref = R0;
    end
end

% check restart. 9-5-17
if ismember('newTend', who('-file', 'parameters.mat'))
    load('parameters.mat', 'newTend');
    if ~isempty(newTend)
        if newTend>tspan(2)
            oldTend = tspan(2);
            tspan(2) = newTend;
            save('parameters.mat', 'tspan', '-append');
            disp(['*** Restarting from T=' num2str(oldTend) ' up to T=' num2str(newTend)]);
            isRestart = true;
        else
            error(['newT=' num2str(newTend) ' < tspan(2)=' num2str(tspan(2))]);
            return
        end
    end
end

% for recording
recNames = {'Y','V','P','VT','YR','B','C','sr2','st2','LA','LAMBDA','TMP','lamBs','PD'};
if exist('isRestart','var') && isRestart
    oldColumns = (oldTend-tspan(1))/dt;
    newColumns = (newTend-oldTend)/dt;
    for i=1:length(recNames)
        switch recNames{i}
            case 'Y'
                load('solution.mat','y'); Y = [y(:,1:oldColumns), zeros(Nr, newColumns)]; clear y
            case 'V'
                load('solution.mat','v'); V = [v(:,1:oldColumns), zeros(Nr, newColumns)]; clear v
            case 'P'
                load('solution.mat','p'); P = [p(:,1:oldColumns), zeros(Nr, newColumns)]; clear p
            case 'sr2'
                load('solution.mat','p'); load('solution.mat','radial');
                sr2 = [radial(:,1:oldColumns)+p(:,1:oldColumns), zeros(Nr, newColumns)]; clear p radial
            case 'st2'
                load('solution.mat','p'); load('solution.mat','hoop');
                st2 = [hoop(:,1:oldColumns)+p(:,1:oldColumns), zeros(Nr, newColumns)]; clear p hoop
            otherwise
                eval(['load(''solution.mat'',''' recNames{i} ''');']);
                eval([recNames{i} ' = [' recNames{i} '(:,1:oldColumns), zeros(Nr, newColumns)];']);
        end
    end
else
    for i=1:length(recNames)
        eval([recNames{i} ' = zeros(Nr, int32(numFrames/record_every+1));']);
    end
end
hoop_warned = false;
disp_progress = 10; % display progress every this many frames. 0 to disable
if exist('radial_grid_search.lock','file'), disp_progress = 0; end
debug_output = 0;

% p, v iteration
use_nonlinear = 0;
max_vp_ite = 100; % max iterations of v,p
tol = 1e-6;
fsolve_ops = optimoptions(@fsolve,'Display','off','algorithm','trust-region');

%% nonlinear solver
    function F = fvp(x)
        % x = [v; p]
        F = ones(size(x));
        
        % v
        lam = get_lambda(sr-x(Nr+1:2*Nr),st-x(Nr+1:2*Nr));
        lamA = get_lambdaA(sr-x(Nr+1:2*Nr),st-x(Nr+1:2*Nr));
        F(1) = x(1); % v(0) = 0
        F(2:Nr) = (2*dx+r(2:Nr)).*x(2:Nr) - r(2:Nr).*x(1:Nr-1) ...
            - dx*R(n-1)*r(2:Nr).*(lam(2:Nr).*c(2:Nr)-lamA(2:Nr));
        
        % p
        F(Nr+1:2*Nr-1) = x(Nr+2:2*Nr)-x(Nr+1:2*Nr-1) ...
            + dx*scale_v*R(n-1)*x(1:Nr-1) - dx*fr(2:Nr);
        F(2*Nr) = x(2*Nr) - pEnd; % p(1) = pEnd
    end

%% feedback functions
    function lambdaA = get_lambdaA(sr2,st2)
        fb = st2;
%         fb = (sr2+2*st2)/3; % include feedback by radial stress
        fs = (abs(fb-s0cA).*(fb<s0cA)).^nlamA; 
        lambdaA = lambdaA_base + lambdaA_A * fcA*fs./(1+fcA*fs);
    end

    function lambda = get_lambda(~,st2)
%         fs = (lambda_A*abs(st2-s0c)).^nlam.*(st2<s0c);
%         lambda = lambda_base ./ (1 + fs);
        
        % Morgan pressure release
        lambda = lambda_base./(1+gLamMns*(abs(st2).*(st2<0)).^nLamMns)...
            +lambda_max*(st2 >0).*(gLamPls*(abs(st2)).^nLamPls)./...
            (1+gLamPls*(abs(st2)).^nLamPls);
%         fs = (abs(st2-s0cA).*(st2>0)).^nlam;
%         lambda = lambda_base + lambda_mr * lambda_A*fs./(1+lambda_A*fs);
    end

    function lamB = get_lamB(~,st2)
        % fs = (abs(st2-s0cL)).^nL.*(st2<s0cL);
        % L = Lbase ./ (1 + fs);
%         fs = gamma_B*(abs(st2-s0cL)).^nL.*(st2<s0cL);
%         lamB = c_lamB ./ (1+fs);
        lamB = 1;%
    end

    function c = get_c(R,r)
%         Lc = Lbase*sqrt(c_lam/(c_lam+c_lamB));
%         c = (c_lam*sinh(R*r/Lc)/sinh(R/Lc)./r + c_lamB) / (c_lam+c_lamB);
%         c(1) = (c_lam*R/(sinh(R/Lc)*Lc) + c_lamB) / (c_lam+c_lamB);
        c = sinh(R*r./Lbase)./(r.*sinh(R./Lbase)); c(1) = R/(sinh(R/Lbase)*Lbase);
    end

%% initialize or restart
r = linspace(0,1,Nr)';
Nt = (tspan(2)-tspan(1))/dt;
dx = r(2)-r(1);
if exist('isRestart','var') && isRestart
    load('solution.mat','R');
    R = [R(1:oldColumns); zeros(newColumns,1)];
    y = Y(:,oldColumns-2);
    nStart = oldColumns-1;
else
    R = zeros(Nt,1); R(1) = R0;
    if ismember('yp',who('-file','parameters.mat'))
        load('parameters.mat','yp');
        y = (r+r.^yp)/2*R0;
        disp(['*** pre-stress p=' num2str(yp)]);
    else
        y = R0*r;
    end
    Y(:,1) = y;
    
    nStart = 2;
    residual = zeros(Nt,1); residual(1) = NaN;
end

%% finite difference loop
for n=nStart:Nt
    currTime = tspan(1)+dt*(n-1);
    if disp_progress && mod(currTime-tspan(1),disp_progress)==0, disp(['T=' num2str(currTime)]); end
    if (match_fig1b && abs(currTime-match_fig1b)<dt) % Helminger Fig 1b, pressure release
        if cH>0
            if ~exist('radial_grid_search.lock','file')
%                 disp('*** pressure release ***');
            end
            load('parameters.mat','tumorID');
            switch tumorID
                case 0, load('helminger_fig1b.mat','f0');
                case 1, load('helminger_fig1b.mat','f07'); f0 = f07;
                case 2, load('helminger_fig1b.mat','f10'); f0 = f10;
            end
            save('parameters.mat','f0','-append');
            cH = 0; 
        end
        if pBar>0
            disp('*** pressure release Morgan ***');
            pBar = 0;
        end
    end
    
    %% calculate stress from lagged p. Note: y is current
    y_r = gradient(y, dx);
    if ~exist('p','var'), p = 0; end % assuming no feedback at T=0
    if with_G_incompatibility
        sr = cT*(y./(r.*y_r)).^(4/3)+cT*mu*(y./(r.*y_r)).^(-4/3); sr(1) = cT+cT*mu;
        st = cT*(y_r.*r./y).^(2/3)+cT*mu*(y_r.*r./y).^(-2/3); st(1) = cT+cT*mu;
    else
        sr = cT/R(n-1)*(y./r).^2./y_r; sr(1) = cT/R(n-1)*y_r(1);
        st = cT/R(n-1)*y_r;
    end
    
    lambda = get_lambda(sr-p,st-p);
    lambdaA = get_lambdaA(sr-p,st-p);
    lamB = get_lamB(sr-p,st-p);
    c = get_c(R(n-1),r);
    
    %% update v and p together
    if debug_output, tic; end
    gel_stress = -cH/2*(5-R0_ref*(R0_ref^3+4*R(n-1)^3)/(R(n-1)^4)) ...
        - cH2*(1-(2*R(n-1)^3-R0_ref^3)/(R(n-1)^2*R0_ref));%-R(n)*rs(R(n),lambda,lambdaA,L);
    y_rr = 4*del2(y, dx);
    if with_G_incompatibility
        sigma_r3 = -cT/R(n-1)*(4*(y./(r.*y_r)).^(1/3).*(-r.*y_r.^2+y.*(y_r+r.*y_rr))./(3*r.^2.*y_r.^2));% sigma_r3(1) = sigma_r3(2);
%         sigma_r3 = gradient(sr,dx)/R(n-1);
        sigma_r3 = sigma_r3 .* (1-mu./sr.^2);
        pEnd = 2*gamma_/R(n-1) + cT*(y(end)/y_r(end))^(4/3) + cT*mu*(y(end)/y_r(end))^(-4/3) - gel_stress + pBar;
    else
        sigma_r3 = 2*cT./(R(n-1)^2*y_r).*(y.*y_r./(r.^2)-y.^2./(r.^3))-cT./(R(n-1)^2*y_r.^2).*(y./r).^2.*y_rr; sigma_r3(1) = sigma_r3(2);
        pEnd = 2*gamma_/R(n-1) + cT/R(n-1)*y(end)^2/y_r(end) - gel_stress + pBar;
    end
    fr = R(n-1)*sigma_r3 + 2./r.*(sr - st);
    N = 2*Nr; % 2x large matrix
    diagV = [1; 2*dx+r(2:end)]; subDiagV = -r(2:end); % backward difference. forward gives singular matrix.
    diagP = [-ones(Nr-1,1); 1]; supDiagP = ones(Nr-1,1); % forward difference
    VPcoupled = [ones(Nr-1,1)*scale_v*R(n-1)*dx; 0]; % lower left submatrix. the p_r+R*v bit
    diag = [diagV; diagP];
    supDiag = [zeros(Nr,1); supDiagP];
    subDiag = [subDiagV; zeros(Nr,1)];
    A = sparse(1:N,1:N,diag) + sparse(1:N-1,2:N,supDiag,N,N) + sparse(2:N,1:N-1,subDiag,N,N) ...
        + sparse(Nr+1:N,1:Nr,VPcoupled,N,N);
    
    if ~use_nonlinear
        for i_vp=1:max_vp_ite % iterate a few times to eliminate oscillation
            % RHS of v and p
            vr = R(n-1)*(lambda.*c-lambdaA).*r;
            rhsV = [0; dx*vr(2:end)];
            rhsP = [dx*fr(2:Nr); pEnd];
            rhs = [rhsV; rhsP];
            % solution is [v; p]
            sol = A\rhs;
            v = sol(1:Nr);
            p = sol(Nr+1:end);
            
            lambda = get_lambda(sr-p,st-p);
            lambdaA = get_lambdaA(sr-p,st-p);
            lamB = get_lamB(sr-p,st-p);
            c = get_c(R(n-1),r);
            
            res = norm(fvp([v; p]));
            if res<tol
                break
            else
                if i_vp==max_vp_ite, use_nonlinear = 1; end
            end
        end
    else
        x0 = [v; p];
        sol = fsolve(@fvp, x0, fsolve_ops);
        v = sol(1:Nr);
        p = sol((Nr+1):(2*Nr));

        lambda = get_lambda(sr-p,st-p);
        lambdaA = get_lambdaA(sr-p,st-p);
        lamB = get_lamB(sr-p,st-p);
        c = get_c(R(n-1),r);
    end
    
    residual(n) = norm(fvp([v; p]));
    if debug_output
        fprintf('T=%.2f, residual=%e, toc=%f\n',currTime,residual(n),toc);
    end
        
    if currTime>0.
    end
    
    % anatical solution of v, p_r
    va = R(n-1)./r.^2.*cumsimps(r,(lambda.*c-lambdaA).*r.^2);
    pra = -scale_v*R(n-1)*va + R(n-1)*sigma_r3 + 2./r.*(sr - st);
    
    % pD is the pressure from Darcy's law
    pD = zeros(Nr,1);
    pD(end) = 2*gamma_/R(n-1) - gel_stress + pBar;
    pDrhs = -scale_v*R(n-1)*va;
    for i_pD=Nr-1:-1:1
        pD(i_pD) = pD(i_pD+1) - dx*pDrhs(i_pD);
    end
    
%     save('new','r','p','st','sr','y','v','pra');
    
    %% update R
    R(n) = R(n-1) + dt*v(end);
    
    %% solve y
    drdt = v(end);
    vt = (v-r*drdt)/R(n);
    %     b = beta_base * ((1-beta_w)*lambda.*c + beta_w*lambdaA );
    b = beta_base * ones(size(r)); % constant remodeling. 7-25-17
    ynew = zeros(Nr,1);
    ynew(end) = (y(end) + dt*R(n)*b(end))/(1+dt*b(end));
    
    % upwind
    vp = max(vt, 0); vn = min(vt, 0);
    for i=2:Nr-1
        upw = (vp(i)*(y(i)-y(i-1))+vn(i)*(y(i+1)-y(i)))/dx;
        ynew(i) = (y(i) - dt*upw + dt*b(i)*R(n)*r(i))/(1+dt*b(i));
    end
    
%     % Lax-Wendroff. has oscillations when y changes sharply in space
%     for i=2:Nr-1 
%         ynew(i) = (y(i) - dt/dx/2*vt(i)*(y(i+1)-y(i-1)) ...
%             + dt^2/dx^2/2*vt(i)^2*(y(i+1)-2*y(i)+y(i-1)) ...
%             + dt*b(i)*R(n)*r(i))/(1+dt*b(i));
%     end

    y = ynew;

    %% bookkeeping
    if ~hoop_warned && max(abs(sr))>1e4 && ~exist('radial_grid_search.lock','file') % do not display warning when fitting
        warning(['T=' num2str(currTime) ', sigma_rr-p too large, ' num2str(max(abs(sr)))]);
        hoop_warned = true;
%         break
    end
    
    if mod(n,record_every)==0
        nn = n/record_every;
        Y(:,nn) = y;
        YR(:,nn) = y_r;
        P(:,nn) = p;
        V(:,nn) = v;
        sr2(:,nn) = sr;
        st2(:,nn) = st;
        VT(:,nn) = vt;
        C(:,nn) = c;
        B(:,nn) = b;
        LA(:,nn) = lambdaA;
        LAMBDA(:,nn) = lambda;
        TMP(:,nn) = y./r/R(n);
        lamBs(:,nn) = lamB;
        VA(:,nn) = va;
        PRA(:,nn) = pra;
        CA(:,nn) = c;
        PD(:,nn) = pD;
    end
end
radial = sr2 - P;
hoop = st2 - P;
end

function m = mina(a,b)
if abs(a)<abs(b)
    m = a;
else
    m = b;
end
end
