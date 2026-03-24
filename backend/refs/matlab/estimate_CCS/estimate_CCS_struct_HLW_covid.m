%this script estimates by maximum likelihood the system with output-asset
%price relation and Phillips curve WITHOUT including r_t in the estimation.

% We treat Covid exactly as in HLW to show the comparison.


%% HOUSEKEEPING
 
clc, clear all, close all

warning('off','MATLAB:dispatcher:nameConflict')
path = '/Users/tomyc/Dropbox (MIT)/FCI_targeting/FCIstar/code/empirics';
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
backward_inf = 1; %use the 1 year past inflation as proxy for future inflation. Same as LW
set_covid_2022 = 0; %1 sets the covid dummy to 0 after 2022.



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
    enddate = find(date == 2024.75);
    startdate_HLW_str = '01-Apr-1990';
enddate_HLW_str = '01-Oct-2024';
elseif recent_data == 4
    startdate = find(date == 1961);
    enddate = find(date == 2024.75);
    startdate_HLW_str = '01-Jan-1961';
    enddate_HLW_str = '01-Oct-2024';
end

[nat_output,output_gap_hp] = hpfilter(gdp_tot,Smoothing = 36000); % separate in "natural" and "gap"
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
% param.c_pi = theta(18);


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
  1 0 0 1]'; % rho_pi rho_delta eta c_g
else
   fix_params_step3 = [1  1   ... %a_{y,1}, a_{y,2}, a_{f}
  1    1, ... % b_{pi}, b_{y} 
  0   0 0 1 0, ... % sigma_yt sigma_pi  sigma_ys sigma_pi_t sigma_delta
  0    0    0   0, ... % phi kappa_vec(1,2,3)
  1 0 0 1]'; % rho_pi rho_delta eta c_g 
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

fixed_value_step3 = [0 0   ... %a_{y,1}, a_{y,2}, a_{f}
  0.689    0.08, ... % b_{pi}, b_{y} 
  0   0 0 0 0, ... % sigma_yt sigma_pi  sigma_ys sigma_pi_t sigma_x
  0    1    1   1, ... % phi kappa_vec(1,2,3)
  1 0.7 0.9 1]'; % rho_pi rho_delta eta c_g

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


%% Create datetime vecs for plots

Y_date_plot = floor(dates_use);
M_date_plot = 12*(dates_use-Y_date_plot)+3;
D_date_plot = 30;

dates_plot = datetime(Y_date_plot,M_date_plot,D_date_plot);

Y_date_plot_fci_t = floor(dates_fci_t);
M_date_plot_fci_t = 12*(dates_fci_t-Y_date_plot_fci_t)+3;
D_date_plot_fci_t = 30;

dates_plot_fci_t = datetime(Y_date_plot_fci_t,M_date_plot_fci_t,D_date_plot_fci_t );


 %% Compare our dummy treatment with HLW.
cd([path vintage task '/_results'])
% load 
main_results = readtable("fci_star_results.xlsx");

fci_star_base = main_results{:,3};
y_gap_base = main_results{:,4};


close all

figure('Position',[0 100 1000 600])

subplot(2,1,1)
grid on, box on
set(gca,'FontSize',12)
hold on
plot(dates_plot,fci_use(startdate:enddate),'LineWidth',3,'Color','k','LineStyle','--')
hold on
plot(dates_plot,fci_star_base,'LineWidth',3,'Color','b')
hold on
plot(dates_plot,fci_star,'LineWidth',3,'Color','r')
hold on
plot(dates_plot,zeros(length(dates_plot),1),'LineWidth',0.5,'Color','k','LineStyle','--')
legend({'FCI','$FCI^{*}$ - Base', '$FCI^{*}$ - HLW Covid'},'Interpreter','latex','Location','southwest','FontSize',12)
title('FCI and $FCI^{*}$','Interpreter','latex','FontSize',18)

subplot(2,1,2)
grid on, box on
set(gca,'FontSize',12)
hold on
plot(dates_plot,zeros(length(dates_plot),1),'LineWidth',0.5,'Color','k','LineStyle','--')
hold on
plot(dates_plot,y_gap_base,'LineWidth',3,'Color','b')
hold on
plot(dates_plot,y_gap,'LineWidth',3,'Color','r')
hold on
legend({'','Base','HLW Covid'},'Location','southwest','FontSize',12,'NumColumns',3)
title('Estimated Output Gap','Interpreter','latex','FontSize',18)


if save_fig == 1
    print(['fci_star_and_output_gap_struct_' num2str(floor(date(startdate))) '_' num2str(floor(date(enddate))) '_HLW_covid'],'-dpng')
    print(['fci_star_and_output_gap_struct_' num2str(floor(date(startdate))) '_' num2str(floor(date(enddate))) '_HLW_covid'],'-depsc')
end

