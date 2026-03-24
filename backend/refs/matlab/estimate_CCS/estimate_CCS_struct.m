%this script estimates by maximum likelihood the system with output-asset
%price relation and Phillips curve WITHOUT including r_t in the estimation.

% I implement a more structural version:


%% HOUSEKEEPING
 
clc, clear all, close all

warning('off','MATLAB:dispatcher:nameConflict')
path = '/Users/tomyc/MIT Dropbox/Tomas Caravello/FCI_targeting/FCIstar/code/empirics';
vintage = '/26_01_10';
task = '/estimate_CCS';

addpath([path vintage '/_aux_fun'])
addpath([path '/_data'])
addpath([path '/_data/HLW_data'])
cd([path vintage task])
%% Settings
save_fig = 0;
recent_data = 3;
run_step_1 = 0; % 1 if you want to rerun step 1, 0 if you are fine calibrating the lambda_g
backward_inf = 1; %use the 1 year past inflation as proxt for future inflation.

save_results = 0;
set_covid_2022 = 1; %1 sets the covid dummy to 0 after 2022.
restrict_af = 1; % 0 if you want a_f to be unrestricted during estimation.

%% DATA
load data_hlm.mat

date = data_var(:,1);
gdp_tot = 100*data_var(:,2);
pce_core = data_var(:,3);
expected_pce = data_var(:,4);
ffr = data_var(:,5);
covid_index = data_var(:,6);
if set_covid_2022 == 1
    stop_covid_date = find(date==2023);
    covid_index(stop_covid_date:end)=0;
end
fci = data_var(:,7); %chicago FED FCI
fci_2 = data_var(:,8); %new FCI
fci_2_1yr = data_var(:,9); %new FCI, 1yr
gdp_pot_cbo = 100*data_var(:,10);
y_gap_cbo = gdp_tot-gdp_pot_cbo;
r_rate = ffr -expected_pce;

% also load covid dummies
covid_dummies = covid_dummies_data(:,2:end);

% import FCI-T for plotting purposes

load fci_t_ols_no_smooth.mat
dates_fci_t = dates_fcst;
fci_t_no_smooth = data_fci_target_retrend_OLS;

clear dates_fcst data_fci_target_retrend_OLS

load fci_t_ols_smooth.mat
fci_t_smooth = data_fci_target_retrend_OLS;
clear data_fci_target_retrend_OLS

%% Settings

if recent_data == 1
startdate = find(date == 1990);
enddate = find(date == 2019.75);
startdate_HLW_str = '01-Jan-1990';
enddate_HLW_str = '01-Oct-2019';

elseif recent_data == 0
startdate = find(date == 1961);
enddate = find(date == 2019.75);
startdate_HLW_str = '01-Jan-1961';
enddate_HLW_str = '01-Oct-2019';
elseif recent_data == 2
    startdate = find(date == 1975);
    enddate = find(date == 2024.75);
    startdate_HLW_str = '01-Jan-1975';
enddate_HLW_str = '01-Oct-2024';
elseif recent_data == 3
    startdate = find(date == 1990.25);
    enddate = find(date == 2025.5);
    startdate_HLW_str = '01-Apr-1990';
enddate_HLW_str = '01-Jul-2025';
elseif recent_data == 4
    startdate = find(date == 1961);
    enddate = find(date == 2024.75);
    startdate_HLW_str = '01-Jan-1961';
    enddate_HLW_str = '01-Oct-2024';
end

[nat_output,output_gap_hp] =hpfilter(gdp_tot,Smoothing=36000); % separate in "natural" and "gap"
%start p_t such that p_t - p^* is close to zero.
%p_0 = nat_output(startdate-4+n_lags)+(nat_output(startdate-4+n_lags)-nat_output(startdate-5+n_lags))-2*log(0.02);
p_0 = 0;
% from the FCI data, build the price index

if recent_data == 3
    fci_use = fci_2; % demean FCI?
else
    fci_use = fci;
end
n_lags = 4;

XX_fci_aux = lagmatrix(fci_use,1:4);
XX_fci= XX_fci_aux(startdate-4+n_lags:end,:);
YY_fci = fci_use(startdate-4+n_lags:end);

beta_fci = (XX_fci'*XX_fci)\(XX_fci'*YY_fci);
delta_p = YY_fci-XX_fci*beta_fci;
p_index = NaN(length(data_var),1);
p_index(startdate-4+n_lags:end) = -cumsum(delta_p)+p_0; %put a minus so higher prices means boost to GDP
p_index(startdate-1) = p_0;


% start with the same set of variable as HLM. Note that GDP tot is already
% multiplied by 100

vardata_aux = [gdp_tot pce_core fci_use];

% create matrix of observables and lags. For all steps, y_mat is the same.

y_mat_aux = vardata_aux(:,[1 2]);

% crate big x matrix with this order: 

%[y_{t-1},y_{t-2},p_{t-1},\pi_{t-1},\pi_{t-2:4},1,d_{t},d_{t-1},d_{t-2}]


x_mat_aux = NaN(size(y_mat_aux,1),9);
x_mat_aux(2:end,1) = vardata_aux(1:end-1,1); %y_{t-1}
x_mat_aux(3:end,2) = vardata_aux(1:end-2,1); %y_{t-2}
x_mat_aux(4:end,3) = vardata_aux(1:end-3,1); %y_{t-3}
x_mat_aux(2:end,4) = vardata_aux(1:end-1,3); % FCI^{m}_{t-1}
x_mat_aux(2:end,5) = vardata_aux(1:end-1,2); % \pi_{t-1}
x_mat_aux(3:end,6) = movmean(vardata_aux(1:end-2,2),[2 0],'Endpoints','shrink'); %\pi_{t-2:4}
x_mat_aux(:,7) = 1; %cosntant
% add the covid index and two lags
x_mat_aux(:,8) = covid_index; %d_{t}
x_mat_aux(2:end,9) = covid_index(1:end-1); % d_{t-1}
x_mat_aux(3:end,10) = covid_index(1:end-2); % d_{t-2}
x_mat_aux(4:end,11) = covid_index(1:end-3); % d_{t-3}

%add linear trend to incorporate when constant trends are needed.
x_mat_aux(:,12) = 1:size(y_mat_aux,1);
x_mat_aux(2:end,13) = 1:size(y_mat_aux,1)-1;
x_mat_aux(3:end,14) = 1:size(y_mat_aux,1)-2;


y_mat = y_mat_aux(startdate:enddate,:); %this is the same
% create the x's for each step.
if backward_inf == 1
    x_mat_step1 = x_mat_aux(startdate:enddate,[1 2 5 6 8 9 10]);

    % x_t = [y_{t-1},FCI_{t-1}^{m},\pi_{t-1},\pi_{t-2:4},d_{t},d_{t-1}]
    x_mat_step3 = x_mat_aux(startdate:enddate,[1 4 5 6 8 9]);
    covid_d_pos_step3 = [5 6]; %position of covid index in x_mat_step3
else
    x_mat_step1 = x_mat_aux(startdate:enddate,[1 2 3 5 8 9 10 11]);
    x_mat_step3 = x_mat_aux(startdate:enddate,[1 2 3 4 5 8 9 10 11]);
    covid_d_pos_step3 = [6 7 8 9]; %position of covid index in x_mat_step3

end

output_nat_init =  nat_output(startdate-1:-1:startdate-4);
% start in the same output gaps of HLW
%output_nat_init = gdp_tot(startdate-1:-1:startdate-4)+y_gap_cbo(startdate-1:-1:startdate-4);
pi_init = expected_pce(startdate-1);% intial values of trend inflation.
fci_init = [0,0];

options = optimset('Display','iter','MaxIter',20000,'MaxFunEvals',20000);
%% load HLW most recent estimate

T_HLW_y_gap = readtable('Holston_Laubach_Williams_current_estimates.xlsx','Sheet','HLW Estimates');

HLW_dates = T_HLW_y_gap{:,1};
startdate_HLW = find(HLW_dates==startdate_HLW_str);
enddate_HLW = find(HLW_dates==enddate_HLW_str);
HLW_y_gap = T_HLW_y_gap{startdate_HLW:enddate_HLW,15};
HLW_r_star =  T_HLW_y_gap{startdate_HLW:enddate_HLW,11};

HLW_r_gap = 0.5*((r_rate(startdate:enddate) - HLW_r_star)+(r_rate(startdate-1:enddate-1) - T_HLW_y_gap{startdate_HLW-1:enddate_HLW-1,11}));
%% Step 1
 if run_step_1 == 1

settings = struct;
settings.smoother = 1;
settings.covid_dummies  = covid_dummies(startdate:enddate,:);

if backward_inf == 1
    map_theta_to_mats_step1 = @map_theta_to_mats_HLM_step1_covid;
    settings.detrend_y = 1;
    settings.states_smoother = 1:3;

    xi_0 = zeros(3,1);
    xi_0(1:3) = output_nat_init(1:3);

    P_0 = eye(3);
else
    map_theta_to_mats_step1 = @map_theta_to_mats_CCS_struct_step1;
    settings.detrend_y = 2;
    settings.states_smoother = 1:5;
    xi_0 = zeros(6,1);
    xi_0(1:4) = output_nat_init;
    xi_0(5) = pi_init;
    xi_0(6) = 1;

    P_0 = eye(6);
    P_0(5,5) = 0.25;
    P_0(6,6) = 0;

end

settings.xi_0 = xi_0;
settings.P_0 = P_0;
settings.fixed_param = [];

% define function to minimize
% USE THE SAME @map_theta_to_mats_HLM_step1_covid file since it is exactly
% the same as HLW.

f_step1 = @(theta) -likelihood_kalman_covid(y_mat,x_mat_step1,theta,map_theta_to_mats_step1,settings);

if backward_inf == 1

% a_y1 = theta(1);
% a_y2 = theta(2);
% b_pi = theta(3);
% b_y = theta(4);
% g = theta(5);
% sigma_yt = theta(6);
% sigma_pi = theta(7);
% sigma_ys = theta(8); % signal
% phi = theta(9);
% kappa_vec = [1;theta(10:12)];

    x0_step1 = [0.3   0, ... % a_y1 a_y2
    0.5    0.08    0.63,... %b_pi b_y g
    0.3    0.2    0.4,... % sigma_yt, sigma_pi, sigma_ys
    -0.1318    8.6454    2.4786    3.1604]';% phi, kappa_vec
LB_step1 = [-2; -2;...
    0; 0.05; 0;...
    0.01; 0.01; 0.01;...
    -inf;1;1;1]';
UB_step1 = [2; 2;...
    inf; inf; 1;...
    inf; inf; inf;...
    inf;inf;inf;inf];
else
% param.a_y1 = theta(1);
% param.a_y2 = theta(2);
% param.b_pi = theta(3);
% param.b_y = theta(4);
% param.g = theta(5);
% param.sigma_yt = theta(6);
% param.sigma_pi = theta(7);
% param.sigma_ys = theta(8);
% param.sigma_pi_e = theta(9);
% param.phi = theta(10);
% param.kappa_vec = [1;theta(11:13)];
% param.rho_pi = theta(14);

    x0_step1 = [0.3   0, ... % a_y1 a_y2
    0.5    0.08    0.63,... %b_pi b_y g
    0.3    0.2    0.4 0.2,... % sigma_yt, sigma_pi, sigma_ys sigma_pi_e
    -0.1318    8.6454    2.4786    3.1604,...% phi, kappa_vec
    0.95]'; %rho_pi
LB_step1 = [-2; -2;...
    0; 0.05; 0;...
    0.01; 0.01; 0.01;0.1;...
    -inf;1;1;1;...
    -1];
UB_step1 = [2; 2;...
    inf; inf; 1;...
    inf; inf; inf;inf;...
    inf;inf;inf;inf;...
    1];
end


[theta_opt_aux,fval,exitflag,~] = fminsearchbnd(f_step1,x0_step1,LB_step1,UB_step1,options);
estimates_step1_aux = evaluate_kalman_covid(y_mat,x_mat_step1,theta_opt_aux,map_theta_to_mats_step1,settings);

% Do same procedure as in HLW: they first estimate the full thing, get the
% filtered uncertainty, and rerun everything using that as initial value

settings.P_0 = squeeze(estimates_step1_aux.state_one_ahead_pred(1,:,:));
[theta_opt,fval,exitflag,~] = fminsearchbnd(f_step1,theta_opt_aux,LB_step1,UB_step1,options);

% this is not EXACTLY the same because we are evaluating at slightly
% different parameters. check: use the same parameters, get exactly the same!
if backward_inf == 1
    kappa_vec = [1;theta_opt(10:12)];
    kappa_vec_use = settings.covid_dummies*kappa_vec;
else
    param = unpack_params_struct_step1(theta_opt);
    kappa_vec_use = settings.covid_dummies*param.kappa_vec;
end

estimates_step1 = evaluate_kalman_covid(y_mat,x_mat_step1,theta_opt,map_theta_to_mats_step1,settings);

% Use estimates in step 1 to compute lambda_g
pot_smoothed_step_1 = estimates_step1.state_smoothed(:,1)/100;
smoothed_growth =400*(pot_smoothed_step_1(2:end)-pot_smoothed_step_1(1:end-1));
lambda_g = median_unbiased_estimator_stage1_covid(pot_smoothed_step_1,kappa_vec_use);
 else
    lambda_g = 0.0667;
 end

%% Step 3

if backward_inf == 1
    map_theta_to_mats_step3 = @map_theta_to_mats_CCS_struct_step3_backward_inf;
else
    map_theta_to_mats_step3 = @map_theta_to_mats_CCS_struct_step3;
end

settings_step3 = struct;
settings_step3.detrend_y = 0;
settings_step3.smoother = 0;
settings_step3.states_smoother = 1:11; % for these states, the smoother is computed. Use to avoid constants, which screw the algorithm.


% states order:
% xi_t = [y_{t}^{n},y_{t-1}^{n},y_{t-2}^{n},g_{t},g_{t-1},FCI_{t}^{*},FCI_{t-1}^{*},\delta_{t},\delta_{t-1}]

y_n_pos = 1;
g_pos = 4;
fci_pos = 7;
delta_pos = 9;


xi_0 = zeros(11,1);
xi_0(1:3) = output_nat_init(1:3);
xi_0(4) = nat_output(startdate-1)-nat_output(startdate-2);
xi_0(5) = nat_output(startdate-2)-nat_output(startdate-3);
xi_0(6) = nat_output(startdate-3)-nat_output(startdate-4);
xi_0(7:8) = [0;0]; % approximate zero initial FCI*
xi_0(9:11) = [0;0;0];

settings_step3.xi_0 = xi_0;
settings_step3.P_0 = 0.25*eye(11);
%settings_step3.P_0(1,1) = 1^2; %output gap has 1% SD
%settings_step3.P_0(2,2) = 1^2;
%settings_step3.P_0(3,3) = 1^2;
settings_step3.P_0([1 2 3],[1 2 3]) = eye(3);
%settings_step3.P_0([4 5 6],[4 5 6]) = 0.25*ones(3,3);  %quarterly growth rate has 0.5 SD
%settings_step3.P_0([7 8],[7 8]) = 0.25*ones(2,2);
%settings_step3.P_0([9 10 11],[9 10 11]) = 0.25*ones(3,3);
%settings_step3.P_0(8,8) = 1^2; 
%settings_step3.P_0(9,9) = 1^2; 


settings_step3.fixed_param = [lambda_g];
settings_step3.covid_dummies = covid_dummies(startdate:enddate,:);
% define function to minimize
f_step3 = @(theta) -likelihood_kalman_covid(y_mat,x_mat_step3,theta,map_theta_to_mats_step3,settings_step3);


x0_step3 = [0.3   0  ... %a_{y,1}, a_{y,2}
 0.5    0.1, ... % b_{pi}, b_{y} 
 0.4140    0.5011    0.2 0.2 0.5, ... % sigma_yt sigma_pi  sigma_ys sigma_pi_e sigma_delta
 -0.1       7    3    2, ... % phi kappa_vec(1,2,3)
 0.9 0.7 0.8 1]'; % rho_pi rho_delta eta  c_pi

% Fix the parameter of dynamics, so we run the structural model directly.
if set_covid_2022 == 1
fix_params_step3 = [1  1   ... %a_{y,1}, a_{y,2}, a_{f}
  1    1, ... % b_{pi}, b_{y} 
  0   0 0 1 0, ... % sigma_yt sigma_pi  sigma_ys sigma_pi_t sigma_delta
  0    0    1   1, ... % phi kappa_vec(1,2,3)
  1 0 0 restrict_af]'; % rho_pi rho_delta eta alpha
else
   fix_params_step3 = [1  1   ... %a_{y,1}, a_{y,2}, a_{f}
  1    1, ... % b_{pi}, b_{y} 
  0   0 0 1 0, ... % sigma_yt sigma_pi  sigma_ys sigma_pi_t sigma_delta
  0    0    0   0, ... % phi kappa_vec(1,2,3)
  1 0 0 restrict_af]'; % rho_pi rho_delta eta alpha
end

LB_step3 = [-2  -2   ... %a_{y,1}, a_{y,2}, a_{f}
  0    0.025, ... % b_{pi}, b_{y} 
  0.01    0.00001 0.0001 0.00001 0.001, ... % sigma_yt sigma_pi  sigma_ys sigma_pi_t sigma_delta
  -inf    1    1    1, ... % phi kappa_vec(1,2,3)
  0 0 0 0]'; % rho_pi rho_delta eta c_g

UB_step3 = [2  2   ... %a_{y,1}, a_{y,2}, a_{f}
  1    inf, ... % b_{pi}, b_{y} 
  inf   inf inf inf inf, ... % sigma_yt sigma_pi  sigma_ys sigma_pi_t sigma_x
  inf    inf    inf   inf, ... % phi kappa_vec(1,2,3)
  1 1 1 5]'; % rho_pi rho_delta eta c_g
 

% choose here at what value you want to fix them.

fixed_value_step3 = [0 0   ... %a_{y,1}, a_{y,2},
  0.689    0.08, ... % b_{pi}, b_{y} 
  0   0 0 0 0, ... % sigma_yt sigma_pi  sigma_ys sigma_pi_t sigma_x
  0    1    1   1, ... % phi kappa_vec(1,2,3)
  1 0.7 0.9 1]'; % rho_pi rho_delta eta alpha

% fix parameters by imposing the constraint.
for i_p = 1:length(x0_step3)
    if fix_params_step3(i_p) == 1
    LB_step3(i_p) = fixed_value_step3(i_p);
    UB_step3(i_p) = fixed_value_step3(i_p);
    end
end

[theta_opt3_aux,~,~,~] = fminsearchbnd(f_step3,x0_step3,LB_step3,UB_step3,options);
estimates_step3_aux = evaluate_kalman_covid(y_mat,x_mat_step3,theta_opt3_aux,map_theta_to_mats_step3,settings_step3);

settings_step3.P_0 = squeeze(estimates_step3_aux.state_contemporaneous_pred(1,:,:));
%f_step3 = @(theta) -likelihood_kalman_covid(y_mat,x_mat_step3,theta,map_theta_to_mats_step3,settings_step3);
[theta_opt3,fval,exitflag,estim_output_3] = fminsearchbnd(f_step3,theta_opt3_aux,LB_step3,UB_step3,options);
%%
theta_opt_use = theta_opt3;
settings_aux = settings_step3;
settings_aux.smoother = 1;
estimates_step3 = evaluate_kalman_covid(y_mat,x_mat_step3,theta_opt_use,map_theta_to_mats_step3,settings_aux);
loglik_aux = sum(estimates_step3.loglik_vec);
%% Compute parameter SEs and states SEs

compute_se = 1;
if compute_se == 1
N_draws = 5000;
fun_se = @(theta) evaluate_kalman_covid(y_mat,x_mat_step3,theta,map_theta_to_mats_step3,settings_aux);
constraints.UB = UB_step3;
constraints.LB = LB_step3;
%a function that is positive if the constraint is NOT satisfied. This is on
%top of the upper and lower bounds.
constraints.const_f = @(theta) (theta(1)+theta(2)>=1)+(theta(16)>=1) + (theta(15)>=1) + sum(theta(5:9)<0);
se_struct = compute_se_kalman(theta_opt3,constraints,fun_se,N_draws);

% using that, build estimates of the states_se
se_report = se_struct.se_vec;
state_ses = se_struct.state_ses;
% natural output SE

y_star_se_last = state_ses(end,y_n_pos);
g_se_last = 4*state_ses(end,g_pos); %multiply by 4 to annualize

y_star_se_avg = mean(state_ses(:,y_n_pos));
g_se_avg = 4*mean(state_ses(:,g_pos)); %multiply by 4 to annualize

fci_se_ts = NaN(length(state_ses(:,y_n_pos)),y_n_pos);
for tt=1:length(state_ses(:,y_n_pos))
    P_tt = squeeze(se_struct.state_unc_mat(tt,:,:));
    fci_se_ts(tt) = sqrt(P_tt(fci_pos,fci_pos));
end
fci_se_last = fci_se_ts(end);
fci_se_avg = mean(fci_se_ts);
end
%% Report Results

% param.a_y1 = theta(1);
% param.a_y2 = theta(2);
% param.a_f = theta(3);
% param.b_pi = theta(4);
% param.b_y = theta(5);
% param.sigma_yt = theta(6);
% param.sigma_pi = theta(7);
% param.sigma_ys = theta(8);
% param.sigma_pi_e = theta(9);
% param.sigma_delta = theta(10);
% param.phi = theta(11);
% param.kappa_vec = [1;theta(12:14)];
% param.rho_pi = theta(15);
% param.rho_delta = theta(16);
% param.eta = theta(17);
param = unpack_params_struct_step3(theta_opt3);


% Define the parameters
parameter_names = {'$a_{y,1}$', '$a_{y,2}$', '$b_{\pi}$', '$b_{y}$', '$\sigma_{\tilde{y}}$', ...
              '$\sigma_{\pi}$', '$\sigma_{y^{*}}$','$\sigma_{\pi^{e}}$','$\sigma_{\delta}$',...
                '$\phi$','$\kappa_{2020}$','$\kappa_{2021}$', '$\kappa_{2022}$',...
              '$\rho_{\pi}$','$\rho_{\delta}$','$\eta$','$c_{g}$','$a(\eta)$'};
% Table 1
% compute SE of a(eta) using delta method.

a_prime = -(param.a_f^2)*(1+2*param.eta+3*param.eta^2);
if compute_se == 1
a_se = abs(a_prime*se_report(16));
else
a_se = NaN(1,1);    
end
order = 1:18;
% clc
% % Print the LaTeX table
% fprintf('\\begin{table}[h!]\n');
% fprintf('\\centering\n');
% fprintf('\\begin{tabular}{|l|cc|}\n');
% fprintf('\\hline\n');
% fprintf('Parameter & Estimate & Standard Error \\\\\n');
% fprintf('\\hline\n');
% 
% for i = 1:length(parameters)
%     j = order(i); 
%     fprintf('%s & %5.3f & %5.3f \\\\\n', parameters{j}, theta_opt3(j), se_report(j));
% end
% fprintf('\\hline\n');
% fprintf('\\end{tabular}\n');
% fprintf('\\caption{Estimated Parameters and standard errors. Log-likelihood = %5.3f}\n',fval);
% fprintf('\\label{tab:parameters}\n');
% fprintf('\\end{table}\n');
% 
% % Table 2
clc
fprintf('\\begin{table}[h!]\n');
fprintf('\\centering\n');
fprintf('\\begin{tabular}{|l|cc|}\n');
fprintf('\\hline\n');
fprintf('Parameter & Estimate & Standard Error \\\\\n');
fprintf('\\hline\n');
fprintf('$\\lambda_{g}$ & %5.3f & \\\\\n', lambda_g);
for i = 1:length(parameter_names)
    j = order(i); 
    if strcmp(parameter_names{j},'$a(\eta)$')
        fprintf('%s & %5.3f &   %5.3f\\\\\n', parameter_names{j}, -param.a_f,a_se);
    else
    if compute_se == 1
    fprintf('%s & %5.3f & %5.3f \\\\\n', parameter_names{j}, theta_opt3(j), se_report(j));
    else
        fprintf('%s & %5.3f & %5.3f \\\\\n', parameter_names{j}, theta_opt3(j), NaN);
    end
    end
end
fprintf('\\hline \n');
fprintf(' Log-likelihood &  %5.3f & \\\\ \n',-fval);
fprintf('\\hline \n');
if compute_se == 1
fprintf(' S.E (sample avg.) &  & \\\\ \n');
fprintf(' $FCI^{*}$ & %5.3f & \\\\ \n',fci_se_avg );
fprintf(' g (annualized) & %5.3f & \\\\ \n',g_se_avg );
fprintf(' $y^{*}$ & %5.3f & \\\\ \n',y_star_se_avg);
fprintf(' &  & \\\\ \n');
fprintf(' S.E (Final Obs.) &  & \\\\ \n');
fprintf(' $FCI^{*}$ & %5.3f & \\\\ \n',fci_se_last );
fprintf(' g (annualized) & %5.3f & \\\\ \n',g_se_last );
fprintf(' $y^{*}$ & %5.3f & \\\\ \n',y_star_se_last);
fprintf('\\hline \n');
end
fprintf('\\end{tabular}\n');
fprintf('\\caption{Estimated Parameters and standard errors. }\n');
fprintf('\\label{tab:parameters}\n');
fprintf('\\end{table}\n');


%% Plot Results



% levels
fci_star_true =  estimates_step3.state_contemporaneous(:,fci_pos);
% measured FCI star
fci_star = fci_star_true;
y_natural = estimates_step3.state_contemporaneous(:,y_n_pos)+param.phi*x_mat_step3(:,covid_d_pos_step3(1));
y_natural_smoothed = estimates_step3.state_smoothed(:,y_n_pos)+param.phi*x_mat_step3(:,covid_d_pos_step3(1));
y_star= estimates_step3.state_contemporaneous(:,y_n_pos+1)+estimates_step3.state_contemporaneous(:,g_pos+1)...
    + estimates_step3.state_contemporaneous(:,delta_pos)-param.rho_delta*estimates_step3.state_contemporaneous(:,delta_pos+1)...
    +param.phi*x_mat_step3(:,covid_d_pos_step3(1)); %including covid
y_gap = gdp_tot(startdate:enddate)-y_star;
y_star_smoothed = estimates_step3.state_smoothed(:,y_n_pos+1)+estimates_step3.state_smoothed(:,g_pos+1)...
    + estimates_step3.state_smoothed(:,delta_pos)-param.rho_delta*estimates_step3.state_smoothed(:,delta_pos+1)...
    +param.phi*x_mat_step3(:,covid_d_pos_step3(1)); %including covid
y_gap_smoothed = gdp_tot(startdate:enddate) -y_star_smoothed;
% y_gap_smoothed = gdp_tot(startdate:enddate)-y_star_smoothed;
fci_gap = fci_use(startdate:enddate)-fci_star;

% compute y* SEs
%composition of y* as a function of the states
star_vec = zeros(11,1);
star_vec(2)=1;star_vec(5)=1;star_vec(9)=1;star_vec(10)=-param.rho_delta;
y_star_se_contemp = zeros(length(y_gap),1);

for tt=1:length(y_gap)
    y_star_se_contemp(tt) = sqrt(star_vec'*squeeze(estimates_step3.state_contemporaneous_pred(tt,:,:))*star_vec);
end

%inf_trend = (1-param.rho_pi)*param.c_pi+ param.rho_pi*estimates_step3.state_contemporaneous(:,11);
%inf_trend_smooth =  (1-param.rho_pi)*2+ param.rho_pi*estimates_step3.state_smoothed(:,11);
dates_use = date(startdate:enddate);

%% Compute variance decomposition for FCI*

T_use = 500; % time horizon to which we compute

% Write the law of motion of the states as a VAR(2)
% order [y_{t}^{n}, g_{t},delta_{t},FCI^{*}_{t}]
% Will write in companion form 

% xi_t = AA * xi_{t-1} + CC * epsilon

% matrix of VAR coefficients

AA = zeros(8,8); %first 4x4 is matrix on first lag, second part for second lag

AA(1,1) = 1; AA(1,2) = 1;
AA(2,2) = 1;
AA(3,3) = param.rho_delta;
% FCI* equation
AA(4,1) = (1/param.a_f) *(-param.eta); AA(4,5) = (1/param.a_f) * param.eta;
AA(4,2) = (1/param.a_f) * (1-param.eta)*param.c_g; AA(4,6) =  (1/param.a_f) * param.eta*param.c_g;
AA(4,3) = (1/param.a_f)*(param.eta+(1-param.rho_delta)*param.rho_delta); AA(4,7) = (1/param.a_f)*(-param.eta*param.rho_delta); 
AA(4,4) = param.eta;
AA(5:8,1:4) = eye(4);

CC = zeros(8,8);
CC(1,1) = param.sigma_ys;
CC(2,2) = param.sigma_ys*lambda_g*param.c_g;
CC(3,3) = param.sigma_delta;
CC(4,1) = (1/param.a_f)*param.sigma_ys; CC(4,2) = (1/param.a_f)*param.sigma_ys*lambda_g*param.c_g; CC(4,3) =  param.sigma_delta*(1/param.a_f)*(-(param.eta + param.rho_delta));

% compute IRFs with respect to the three shocks.

IRF_states = zeros(8,T_use,3);

for n_shock = 1:3
    IRF_states(:,1,n_shock) = CC(:,n_shock);
    for t=2:T_use
        IRF_states(:,t,n_shock) = AA*squeeze(IRF_states(:,t-1,n_shock));
    end
end

IRF_y_star = squeeze(IRF_states(4,:,1))';
IRF_g = squeeze(IRF_states(4,:,2))';
IRF_delta = squeeze(IRF_states(4,:,3))';

total_var = cumsum(IRF_y_star.^2 + IRF_g.^2 + IRF_delta.^2);

share_delta = cumsum(IRF_delta.^2)./total_var;
share_y_star = cumsum(IRF_y_star.^2)./total_var;
share_g = 1-share_delta-share_y_star;

%%

matrices = map_theta_to_mats_step3(theta_opt3,settings_step3.fixed_param,settings_step3.covid_dummies);

% fitted values
y_fitted = x_mat_step3*(matrices.A) + estimates_step3.state_contemporaneous*(matrices.H);

y_gap_fitted = y_fitted(:,1)-y_star;
inf_fit = y_fitted(:,2);
% fitted inflation

error_inf = y_mat(:,2)-inf_fit;
%%
% Sim
y_gap_sim = y_gap(1)+param.a_f*cumsum(fci_gap);
y_gap_sim_hlw = zeros(length(y_gap_sim),1);
y_gap_sim_hlw(1:2) = HLW_y_gap(1:2);

a_y1_HLW = 1.41652633916052;
a_y2_HLW = -0.483491293620491;
a_r_HLW = -0.0683002396959628;

for tt=3:length(y_gap_sim)
    y_gap_sim_hlw(tt) = a_y1_HLW*y_gap_sim_hlw(tt-1) + a_y2_HLW*y_gap_sim_hlw(tt-2) + a_r_HLW*HLW_r_gap(tt-1); 
end

our_diff = y_gap-y_gap_sim;
HLW_dif = HLW_y_gap-y_gap_sim_hlw;

mat_compare_HLW = [HLW_y_gap y_gap_sim_hlw];

%% Compute historical decomposition

% from the inferred values for the states, we can use the law of motion of
% FCI* to visualize how the observed path of FCI* is explained by each of
% the shocks.

y_n_2 = estimates_step3.state_smoothed(:,y_n_pos);
g = estimates_step3.state_smoothed(:,g_pos);
g_l1 = estimates_step3.state_smoothed(:,g_pos+1);
delta = estimates_step3.state_smoothed(:,delta_pos);
delta_l1 = estimates_step3.state_smoothed(:,delta_pos+1);
delta_l2 = estimates_step3.state_smoothed(:,delta_pos+2);
eps_y_n = estimates_step3.state_smoothed(:,y_n_pos)-(estimates_step3.state_smoothed(:,y_n_pos+1)+estimates_step3.state_smoothed(:,g_pos+1));
eps_y_n_l1 = estimates_step3.state_smoothed(:,y_n_pos+1)-(estimates_step3.state_smoothed(:,y_n_pos+2)+estimates_step3.state_smoothed(:,g_pos+2));
eps_g = estimates_step3.state_smoothed(:,g_pos)-estimates_step3.state_smoothed(:,g_pos+1);
eps_delta = estimates_step3.state_smoothed(:,delta_pos)-param.rho_delta*estimates_step3.state_smoothed(:,delta_pos+1);
eps_delta_l1 = estimates_step3.state_smoothed(:,delta_pos+1)-param.rho_delta*estimates_step3.state_smoothed(:,delta_pos+2);

fci_star_smoothed = estimates_step3.state_smoothed(:,fci_pos);
fci_star_y_n = zeros(length(g),1);
fci_star_g = zeros(length(g),1);
fci_star_delta = zeros(length(g),1);

for tt=2:length(g)
    fci_star_y_n(tt) = param.eta*fci_star_y_n(tt-1)+(1/param.a_f)*(eps_y_n(tt)-param.eta*eps_y_n_l1(tt));
    %fci_star_y_n(tt) = param.eta*fci_star_y_n(tt-1)+(1/param.a_f)*(y_n_2(tt)-param.eta*eps_y_n_l1(tt));
    fci_star_g(tt) = param.eta*fci_star_g(tt-1)+(1/param.a_f)*(eps_g(tt)+(1-param.eta)*g_l1(tt));
    fci_star_delta(tt) = param.eta*fci_star_delta(tt-1) + (1/param.a_f)*(-eps_delta(tt)*(param.eta + param.rho_delta) ...
        +(param.eta + param.rho_delta*(1-param.rho_delta))*delta_l1(tt)-param.eta*param.rho_delta*delta_l2(tt));
end

init_contrib_delta = (1/param.a_f)*(-eps_delta(1)*(param.eta + param.rho_delta)+...
    (param.eta + param.rho_delta*(1-param.rho_delta))*delta_l1(1)-param.eta*param.rho_delta*delta_l2(1));

init_val = fci_star_smoothed(1)*(param.eta.^((0:length(g)-1)'));
fci_star_hd_error = fci_star_smoothed-fci_star_y_n-fci_star_g-fci_star_delta-init_val;
%plot(fci_star_hd_error) %the difference is the intiial value.. now corrected.

% also compute the dynamics starting from a given initial state, will need
% to assign composition to the initial value.

IRF_growth_init = zeros(8,T_use);
IRF_delta_init = zeros(8,T_use);
IRF_growth_init(2,1) = g_l1(1);
IRF_growth_init(6,1) = 0; %put 0 here. This comes from writing eps_{t}^{n} = y^{n}_{t}-(y^{n}_{t}+g_{t-1}) and it is already accounted for.
IRF_delta_init(3,1) = delta_l1(1);
IRF_delta_init(7,1) = delta_l2(1);
    for t=2:T_use
        IRF_growth_init(:,t) = AA*IRF_growth_init(:,t-1);
        IRF_delta_init(:,t) = AA*IRF_delta_init(:,t-1);
    end

init_val_y_n = (eps_y_n(1)/param.sigma_ys)*IRF_y_star(1:length(init_val));
init_val_g = g(1)/param.a_f;
init_val_delta = init_val-init_val_g-init_val_y_n;

fci_star_long_run = param.c_g*(g/param.a_f);
%% Check if fci_gap predicts FCI

n_lags = 2;

XX_fci_aux = lagmatrix(fci_use,1:n_lags);
YY_fci = fci_use(startdate+n_lags:end);
XX_fci= [ones(length(YY_fci),1),XX_fci_aux(startdate+n_lags:end,:)];

beta_fci = (XX_fci'*XX_fci)\(XX_fci'*YY_fci);
u_reg = YY_fci-XX_fci*beta_fci;
var_u_reg = u_reg'*u_reg/(length(u_reg)-length(beta_fci));
V_mat_reg = var_u_reg*inv(XX_fci'*XX_fci);
se_reg = sqrt(diag(V_mat_reg));
t_stats_reg = beta_fci./se_reg;

% het robust SEs
V_mat_het_robust = inv(XX_fci'*XX_fci)*(XX_fci'*((u_reg*u_reg').*eye(length(u_reg)))*XX_fci)*inv(XX_fci'*XX_fci);
se_reg_robust = sqrt(diag(V_mat_het_robust));
t_stats_robust = beta_fci./se_reg_robust;

%essentially clean errors.
errors_acf = autocorr(u_reg);
%%
% add lagged FCI_gap as predictor.
XX_fci_2 = [XX_fci fci_gap(n_lags:end-1) fci_gap(n_lags-1:end-2)];
beta_fci_2 = (XX_fci_2'*XX_fci_2)\(XX_fci_2'*YY_fci);
u_reg_2 = YY_fci-XX_fci_2*beta_fci_2;
V_mat_het_robust_2 = inv(XX_fci_2'*XX_fci_2)*(XX_fci_2'*((u_reg_2*u_reg_2').*eye(length(u_reg_2)))*XX_fci_2)*inv(XX_fci_2'*XX_fci_2);
se_reg_robust_2 = sqrt(diag(V_mat_het_robust_2));

% FCI gap does NOT have additional predictive power beyond what is
% contained in FCI levels.
t_stats_robust_2 = beta_fci_2./se_reg_robust_2;

%% Compute implied effects on output of FCI gaps and r gaps.

HLW_a_r = -0.0683002396959628;

%they actually have the moving average of two periods, but proxy with this.
HLW_gap_on_ouput = -HLW_a_r*HLW_r_gap;
fci_gap_on_output = -param.a_f*fci_gap;

%% Create datetime vecs for plots

Y_date_plot = floor(dates_use);
M_date_plot = 12*(dates_use-Y_date_plot)+3;
D_date_plot = 30;

dates_plot = datetime(Y_date_plot,M_date_plot,D_date_plot);

Y_date_plot_fci_t = floor(dates_fci_t);
M_date_plot_fci_t = 12*(dates_fci_t-Y_date_plot_fci_t)+3;
D_date_plot_fci_t = 30;

dates_plot_fci_t = datetime(Y_date_plot_fci_t,M_date_plot_fci_t,D_date_plot_fci_t );
 %% Replication of HLW Figure 1: very good!
cd([path vintage task '/_results'])
close all

figure('Position',[0 100 1000 600])

subplot(2,1,1)
grid on, box on
set(gca,'FontSize',12)
hold on
plot(dates_plot,fci_use(startdate:enddate),'LineWidth',3,'Color','k','LineStyle','--')
hold on
plot(dates_plot,fci_star,'LineWidth',3,'Color','b')
hold on
plot(dates_plot,zeros(length(dates_plot),1),'LineWidth',0.5,'Color','k','LineStyle','--')
legend({'FCI','$FCI^{*}$'},'Interpreter','latex','Location','southwest','FontSize',12)
title('FCI and $FCI^{*}$','Interpreter','latex','FontSize',18)

subplot(2,1,2)
grid on, box on
set(gca,'FontSize',12)
hold on
plot(dates_plot,zeros(length(dates_plot),1),'LineWidth',0.5,'Color','k','LineStyle','--')
hold on
plot(dates_plot,y_gap,'LineWidth',3,'Color','b')
hold on
%plot(dates_plot,y_gap_smoothed,'LineWidth',3,'Color','b','LineStyle','--')
%hold on
plot(dates_plot,HLW_y_gap,'LineWidth',3,'Color','r','LineStyle','--')
hold on
plot(dates_plot,y_gap_cbo(startdate:enddate),'LineWidth',3,'Color','k','LineStyle','--')
hold on
legend({'','Model','HLW','CBO'},'Location','southwest','FontSize',12,'NumColumns',3)
title('Estimated Output Gap','Interpreter','latex','FontSize',18)


if save_fig == 1
    print(['fci_star_and_output_gap_struct_' num2str(floor(date(startdate))) '_' num2str(floor(date(enddate)))],'-dpng')
    print(['fci_star_and_output_gap_struct_' num2str(floor(date(startdate))) '_' num2str(floor(date(enddate)))],'-depsc')
end
%%

f = figure(2);
f.Position = [0 100 1000 600];
subplot(2,1,1)
grid on, box on
% output gaps
set(gca,'FontSize',12)
hold on
plot(dates_plot,y_gap(:,1),'Color','k','LineWidth',3)
hold on
plot(dates_plot,y_gap_fitted,'r','LineWidth',3)
hold on
plot(dates_plot,zeros(length(dates_plot),1),'LineWidth',0.5,'Color','k','LineStyle','--')
legend({'Estimated','Fitted'},'Location','southwest')
title('Estimated vs Simulated Output Gap','Interpreter','latex','FontSize',16)

% inflation
subplot(2,1,2)
grid on, box on
% output gaps
set(gca,'FontSize',12)
hold on
plot(dates_plot,y_mat(:,2),'Color','k','LineWidth',3)
hold on
plot(dates_plot,inf_fit,'Color','r','LineWidth',3)
hold on
if backward_inf ~= 1
plot(dates_plot,inf_trend,'Color','b','LineWidth',3)
hold on
end
plot(dates_plot,zeros(length(dates_plot),1),'LineWidth',0.5,'Color','k','LineStyle','--')
ylim([-1 9])
legend({'Actual','Fitted'},'Location','northwest','Interpreter','latex')
title('Actual vs Simulated Inflation','Interpreter','latex','FontSize',16)

if save_fig == 1
    print(['actual_vs_fitted_struct_' num2str(floor(date(startdate))) '_' num2str(floor(date(enddate)))],'-dpng')
    print(['actual_vs_fitted_struct_' num2str(floor(date(startdate))) '_' num2str(floor(date(enddate)))],'-depsc')
end

%%

init_gfc = find(dates_use == 2007);

figure(3)
% plot(dates_use(init_gfc:end),y_mat(init_gfc:end,1)-y_mat(init_gfc,1),'Color','k','LineWidth',3)
% hold on
% plot(dates_use(init_gfc:end),y_star(init_gfc:end)-y_mat(init_gfc,1),'b','LineWidth',3)
% hold on
grid on, box on
% output gaps
set(gca,'FontSize',12)
hold on
plot(dates_plot(init_gfc:end),y_mat(init_gfc:end,1)/100,'Color','k','LineWidth',3)
hold on
plot(dates_plot(init_gfc:end),y_natural(init_gfc:end)/100,'b','LineWidth',3)
legend({'Actual','Potential'},'Location','southeast')
ylabel('log(GDP)','Interpreter','latex')
title('Actual vs Potential Output','Interpreter','latex','FontSize',16)
ylim([9.6 10.2])
if save_fig == 1
    print(['gdp_vs_natural_struct_' num2str(floor(date(startdate))) '_' num2str(floor(date(enddate)))],'-dpng')
    print(['gdp_vs_natural_struct_' num2str(floor(date(startdate))) '_' num2str(floor(date(enddate)))],'-depsc')
end

%%

f_4 = figure();
f_4.Position = [0 100 1000 600];

grid on, box on
% output gaps
set(gca,'FontSize',14)
hold on
plot(dates_plot,fci_use(startdate:enddate),'LineWidth',3,'Color','k')
hold on
plot(dates_plot,fci_star,'LineWidth',3,'Color','b')
hold on
%plot(dates_plot,fci_star_smoothed,'LineWidth',3,'Color','b','LineStyle',':')
%hold on
plot(dates_plot_fci_t,fci_t_smooth,'LineWidth',3,'Color','r')
hold on
%plot(dates_plot,fci_zero_y_gap,'LineWidth',3,'Color','r')
plot(dates_plot,zeros(length(dates_plot),1),'LineWidth',0.5,'Color','k','LineStyle','--')
legend({'FCI','$FCI^{*}$','$\overline{FCI}$',''},'Interpreter','latex','Location','southwest')
title('Financial Conditions Index: Data (Dashed), Natural Level (Blue) and Target (red)','Interpreter','latex','FontSize',18)
if save_fig == 1
    print(['fci_star_struct_' num2str(floor(date(startdate))) '_' num2str(floor(date(enddate)))],'-dpng')
    print(['fci_star_struct_' num2str(floor(date(startdate))) '_' num2str(floor(date(enddate)))],'-depsc')
end

%% Save this data in a excel file
fci_t = NaN(length(fci_star),1);
fci_t(1:length(fci_t_smooth)) = fci_t_smooth;
TT = table(dates_plot,fci_use(startdate:enddate),fci_star,y_gap,fci_t);
TT_2 = renamevars(TT,["dates_plot","Var2"],["Date","FCI"]);


if save_results == 1
filename = 'fci_star_results.xlsx';
writetable(TT_2,filename,'Sheet',1)
end

%% Plot FCI* vs HLW r*
f_5 = figure();
f_5.Position = [0 100 1000 600];

grid on, box on
% output gaps
set(gca,'FontSize',14)
hold on
plot(dates_plot,(fci_star-mean(fci_star))/std(fci_star),'LineWidth',3,'Color','b')
hold on
plot(dates_plot,(HLW_r_star-mean(HLW_r_star))/std(HLW_r_star),'LineWidth',3,'Color','r')
hold on
%plot(dates_plot,fci_zero_y_gap,'LineWidth',3,'Color','r')
plot(dates_plot,zeros(length(dates_plot),1),'LineWidth',0.5,'Color','k','LineStyle','--')
hold on
ylabel('Std. Devs from Mean','Interpreter','latex')
legend({'$FCI^{*}$','$r^{*}$'},'Interpreter','latex','Location','southwest')
title('Standardized $FCI^{*}$ and $r^{*}$','Interpreter','latex','FontSize',18)

if save_fig == 1
    print(['r_star_vsfci_star_struct_' num2str(floor(date(startdate))) '_' num2str(floor(date(enddate)))],'-dpng')
    print(['r_star_vsfci_star_struct_' num2str(floor(date(startdate))) '_' num2str(floor(date(enddate)))],'-depsc')
end

%% Plot: current effect on output of gaps

f_6 = figure();
f_6.Position = [0 100 1000 600];
grid on, box on
% output gaps
set(gca,'FontSize',14)
hold on
plot(dates_plot,fci_gap_on_output,'LineWidth',3,'Color','b')
hold on
plot(dates_plot,HLW_gap_on_ouput,'LineWidth',3,'Color','r')
hold on
%plot(dates_plot,fci_zero_y_gap,'LineWidth',3,'Color','r')
plot(dates_plot,zeros(length(dates_plot),1),'LineWidth',0.5,'Color','k','LineStyle','--')
legend({'$FCI$ Gap','$r$ Gap'},'Interpreter','latex','Location','southwest')
title('Impact Effect on Output Gap','Interpreter','latex','FontSize',18)

if save_fig == 1
    print(['effect_on_gaps_struct_' num2str(floor(date(startdate))) '_' num2str(floor(date(enddate)))],'-dpng')
    print(['effect_on_gaps_struct_' num2str(floor(date(startdate))) '_' num2str(floor(date(enddate)))],'-depsc')
end


%% Variance decomposition

colors_use = [0 0.4470 0.7410; 0.85 0.3250 0.0980; 0.4660 0.6740 0.1880; 0.5 0.5 0.5];
alpha_color = 0.75;

YY_var = [share_y_star,share_delta,share_g];
var_hor_plot = 51;
f_7 = figure();
f_7.Position = [0 100 1000 600];
set(gca,'FontSize',14)
hold on
area(0:(var_hor_plot-1),YY_var(1:var_hor_plot,:))
hold on
colororder(colors_use*alpha_color+(1-alpha_color)*ones(size(colors_use)))
legend({'Natural Output','Demand','Growth'},'Interpreter','latex','location','southwest')
xlabel('Quarters')
title('Share of $Var_{t}[FCI_{t+h}^{*}]$ explained by each shock','Interpreter','latex')

if save_fig == 1
    print(['var_decomp_struct_' num2str(floor(date(startdate))) '_' num2str(floor(date(enddate)))],'-dpng')
    print(['var_decomp_struct_' num2str(floor(date(startdate))) '_' num2str(floor(date(enddate)))],'-depsc')
end

%% Composition of FCI*
T_hor_use = 100;
% Assignment of initial values: - the effect of natural output shocks is
% direct. 
% for the split of the remainder, use the contribution to forecast variance
% decomposition between delta and g shocks 25 years in the future (100
% quarters). We do this because the long run variance does not exist.

%delta_init_use = share_delta(T_hor_use)/(share_g(T_hor_use)+share_delta(T_hor_use));
delta_init_use  = 0;
%contrib_mat = [fci_star_y_n fci_star_delta fci_star_g init_val];
%contrib_mat = [(fci_star_y_n + init_val_y_n) (fci_star_delta + delta_init_use*init_val_others) (fci_star_g + (1-delta_init_use)*init_val_others)];
%contrib_mat = [(fci_star_y_n + init_val_y_n) (fci_star_delta + init_val_delta) (fci_star_g + init_val_g)];
contrib_mat = [fci_star_y_n fci_star_delta fci_star_g+init_val];
component_names = {'Natural Output','Demand','Growth','$FCI^{*}$ - Two Sided','$FCI^{*}$ - One Sided'};
f_7 = figure();
f_7.Position = [0 -100 1200 700];
set(gca,'FontSize',14)
hold on
b = bar(dates_use,contrib_mat,"stacked");
hold on
colororder(colors_use*alpha_color+(1-alpha_color)*ones(size(colors_use)))
hold on
plot(dates_use,fci_star_smoothed,'LineWidth',3,'color','black')
hold on
plot(dates_use,fci_star,'LineWidth',3,'color','black','LineStyle','--')
hold on
legend(component_names,'Location','southwest','NumColumns',2,'Interpreter','latex')
title('Historical Decomposition of $FCI^{*}$','Interpreter','latex')
if save_fig == 1
    print(['hist_decomp_struct_' num2str(floor(date(startdate))) '_' num2str(floor(date(enddate)))],'-dpng')
    print(['hist_decomp_struct_' num2str(floor(date(startdate))) '_' num2str(floor(date(enddate)))],'-depsc')
end

%% Plot one and two sided estimates of the states
if compute_se == 1

states_pos = [y_n_pos;g_pos;fci_pos;delta_pos];
states_title = {'$y_{n}$','$g$','$FCI^{*}$','$\delta$'};
f_7 = figure();
f_7.Position = [0 100 1400 800];

bands_conf = 90;
crit_val = norminv(1-(1-bands_conf/100)/2);


for nn=1:4
    subplot(2,2,nn)
grid on, box on
set(gca,'FontSize',12)
hold on
plot(dates_plot,estimates_step3.state_contemporaneous(:,states_pos(nn))-crit_val*se_struct.state_ses_contemp(:,states_pos(nn)),'LineWidth',1,'Color','b','LineStyle','--')
hold on
plot(dates_plot,estimates_step3.state_contemporaneous(:,states_pos(nn))+crit_val*se_struct.state_ses_contemp(:,states_pos(nn)),'LineWidth',1,'Color','b','LineStyle','--')
hold on
plot(dates_plot,estimates_step3.state_contemporaneous(:,states_pos(nn)),'LineWidth',3,'Color','b')
hold on
plot(dates_plot,estimates_step3.state_smoothed(:,states_pos(nn)),'LineWidth',3,'Color','k')
hold on
if nn == 1
    legend({[num2str(bands_conf) ' perc. conf. bands'],'','One-Sided','Two-Sided'},'Location','southeast','FontSize',12)
end

title(states_title{nn},'Interpreter','latex','FontSize',18)
end
if save_fig == 1
    print(['one_sided_two_sided_' num2str(floor(date(startdate))) '_' num2str(floor(date(enddate)))],'-dpng')
    print(['one_sided_two_sided_' num2str(floor(date(startdate))) '_' num2str(floor(date(enddate)))],'-depsc')
end
%%

figure('Position',[0 50 1000 1000])

subplot(2,1,1)
grid on, box on
set(gca,'FontSize',12)
hold on
plot(dates_plot,fci_star,'LineWidth',3,'Color','b')
hold on
plot(dates_plot,fci_star_smoothed,'LineWidth',3,'Color','k')
hold on
hold on
plot(dates_plot,fci_star-crit_val*se_struct.state_ses_contemp(:,fci_pos),'LineWidth',1,'Color','b','LineStyle','--')
hold on
plot(dates_plot,fci_star+crit_val*se_struct.state_ses_contemp(:,fci_pos),'LineWidth',1,'Color','b','LineStyle','--')
hold on
plot(dates_plot,zeros(length(dates_plot),1),'LineWidth',0.5,'Color','k','LineStyle','--')
hold on
%legend({'FCI','$FCI^{*}$'},'Interpreter','latex','Location','southwest','FontSize',12)
title('$FCI^{*}$','Interpreter','latex','FontSize',18)

subplot(2,1,2)
grid on, box on
set(gca,'FontSize',12)
hold on
plot(dates_plot,y_gap-crit_val*y_star_se_contemp,'LineWidth',1,'Color','b','LineStyle','--')
hold on
plot(dates_plot,y_gap+crit_val*y_star_se_contemp,'LineWidth',1,'Color','b','LineStyle','--')
hold on
plot(dates_plot,y_gap,'LineWidth',3,'Color','b')
hold on
plot(dates_plot,y_gap_smoothed,'LineWidth',3,'Color','k')
hold on
plot(dates_plot,zeros(length(dates_plot),1),'LineWidth',0.5,'Color','k','LineStyle','--')
hold on
%legend({'','Model','HLW','CBO'},'Location','southwest','FontSize',12,'NumColumns',3)
title('Output Gap','Interpreter','latex','FontSize',18)
legend({[num2str(bands_conf) ' perc. conf. bands'],'','One-Sided','Two-Sided'},'Location','southwest','FontSize',12)

if save_fig == 1
    print(['fci_output_gap_smoothed_' num2str(floor(date(startdate))) '_' num2str(floor(date(enddate)))],'-dpng')
    print(['fci_output_gap_smoothed_' num2str(floor(date(startdate))) '_' num2str(floor(date(enddate)))],'-depsc')
end
end
