addpath('solver')
addpath('data_as_diameter')
addpath('plotting')

%% parameters for Fig 1
param = [
    1.1;    % lambda_0, eq. 8
    0.2;    % lambda_A0, eq. 9
    40;     % L, eq. 7
    0;      % (not used)
    0;      % (not used)
    0;      % hydrostatic pressure in surrounding material, \bar p
    0;      % stress threshold for feedback to come in, set to 0 in eq. 8
    2;      % m, in eq. 9
    0.6;    % beta, eq. 2
    0;      % (not used)
    0.9;    % Delta_A, eq. 9
    0.1;    % gamma_A, eq. 9
    0;      % c_H, eq. 5
    0;      % tumor ID for experimental data, see radial_set_parameters, line 47
    0;      % (not used)
    0;      % (not used)
    0;      % (not used)
    0;      % gamma_lambda, eq. 8
    2;      % n, eq. 8
    0;      % lambda_max, max effect of positive feedback on lambda (not in paper)
    0;      % gLamPls, positive feedback strength on lambda (not in paper)
    0;      % nLamPls, exponent associated with posotive feedback for lambda (not in paper)
];

%% run simulation
radial_time_evolution(param);

% solution.mat structure:
% 'R' is tumor radius
% solved field functions are Nr*Nt arrays. For example, p(:,1) is the
%   initial condition of pressure, p(end,:) is the boundary solution of
%   pressure for all times.
%   Note: 'radial' and 'hoop' are the total radial/hoop stress.

% solution_rhoc.mat:
% 'RHOC' is the solved rho_c

%% make figures from solution*.mat
radial_plot_single_run;
