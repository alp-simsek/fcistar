%this script estimates by maximum likelihood the system in Holston,
%Laubach, Williams.

% Here I implement their 3-step procedure, that they argue works better for
% econometric reasons. 

%% HOUSEKEEPING
 
clc, clear all, close all

warning('off','MATLAB:dispatcher:nameConflict')
path = '/Users/tomyc/Dropbox (MIT)/FCI_targeting/FCIstar/code/empirics';
vintage = '/25_03_22';
task = '\estimate_HLW';

data_path = '/Users/tomyc/Dropbox (MIT)/FCI_targeting/code/empirics';

addpath([path vintage '/_aux_fun'])
addpath([path '/_data'])
addpath([path '/_data/HLW_data'])
cd([path vintage task])
%% Settings
save_fig = 0;
recent_data = 4;
run_steps_1_2 = 0; %1 if you want to run them, 0 if you want to use the values from HLW directly.
use_HLW_data = 1;% if 1, uses HLW original data. 0 uses my database.

%% DATA
load data_hlm.mat

date = data_var(:,1);
gdp_tot = 100*data_var(:,2);
pce_core = data_var(:,3);
expected_pce = data_var(:,4);
ffr = data_var(:,5);
covid_index = data_var(:,6);

% also load covid dummies
covid_dummies = covid_dummies_data(:,2:end);

% Load HLW original data

TT_HLM = readtable("Holston_Laubach_Williams_current_estimates.xlsx",'Sheet','US input data');
gdp_tot_HLW = nan(length(date),1);
expected_pce_HLW = nan(length(date),1);
pce_core_HLW = nan(length(date),1);
ffr_HLW = nan(length(date),1);
covid_index_HLW = nan(length(date),1);

% start from 5th since this data starts in 1960, not 1959.
gdp_tot_HLW(5:end) = 100*TT_HLM{:,2};
pce_core_HLW(5:end)= TT_HLM{:,3};
expected_pce_HLW(5:end) = TT_HLM{:,4};
ffr_HLW(5:end) = TT_HLM{:,5};
covid_index_HLW(5:end) = TT_HLM{:,6};

% Choose what data to use

if use_HLW_data == 1
    gdp_tot = gdp_tot_HLW;
    pce_core = pce_core_HLW;
    expected_pce = expected_pce_HLW;
    ffr = ffr_HLW;
    covid_index = covid_index_HLW;
end


%% Settings

if recent_data == 1
startdate = find(date == 1990);
enddate = find(date == 2019.75);
elseif recent_data == 0
startdate = find(date == 1961);
enddate = find(date == 2019.75);
elseif recent_data == 2
    startdate = find(date == 1973);
    enddate = find(date == 2024.75);
elseif recent_data == 3
    startdate = find(date == 1991);
    enddate = find(date == 2024.75);
elseif recent_data == 4
    startdate = find(date == 1961);
    enddate = find(date == 2024.75);
end




% start with the same set of variable as HLM. Note that GDP tot is already
% multiplied by 100
vardata_aux = [gdp_tot pce_core (ffr-expected_pce)];

% create matrix of observables and lags. For all steps, y_mat is the same.

y_mat_aux = vardata_aux(:,[1 2]);

% crate big x matrix with this order: 

%[y_{t-1},y_{t-2},r_{t-1},r_{t-2},\pi_{t-1},\pi_{t-2:4},1,d_{t},d_{t-1},d_{t-2}]


x_mat_aux = NaN(size(y_mat_aux,1),10);
x_mat_aux(2:end,1) = vardata_aux(1:end-1,1);
x_mat_aux(3:end,2) = vardata_aux(1:end-2,1);
x_mat_aux(2:end,3) = vardata_aux(1:end-1,3);
x_mat_aux(3:end,4) = vardata_aux(1:end-2,3);
x_mat_aux(2:end,5) = vardata_aux(1:end-1,2);
x_mat_aux(3:end,6) = movmean(vardata_aux(1:end-2,2),[2 0],'Endpoints','shrink');
x_mat_aux(:,7) = 1;
% add the covid index and two lags
x_mat_aux(:,8) = covid_index;
x_mat_aux(2:end,9) = covid_index(1:end-1);
x_mat_aux(3:end,10) = covid_index(1:end-2);
%add linear trend to incorporate when constant trends are needed.
x_mat_aux(:,11) = 1:size(y_mat_aux,1);
x_mat_aux(2:end,12) = 1:size(y_mat_aux,1)-1;
x_mat_aux(3:end,13) = 1:size(y_mat_aux,1)-2;


y_mat = y_mat_aux(startdate:enddate,:); %this is the same
% create the x's for each step.
x_mat_step1 = x_mat_aux(startdate:enddate,[1 2 5 6 8 9 10]);
x_mat_step2 = x_mat_aux(startdate:enddate,1:10);
x_mat_step3 = x_mat_aux(startdate:enddate,[1 2 3 4 5 6 8 9 10]);

% clear this variable because otherwise MATLAB gets confused and thinks
% this is infinity.

% this is used only for starting values. So only use data until 2019Q4 for
% this

%date_precovid = find(date == 2019.75);
[nat_output,output_gap_hp] =hpfilter(gdp_tot,36000); % separate in "natural" and "gap"
if use_HLW_data == 1
    nat_output=[NaN(4,1);nat_output]; %add first 4 NANs so that gdp_tot and nat_output have same length.
end


options = optimset('Display','iter','MaxIter',10000,'MaxFunEvals',10000);
%% Step 1

if run_steps_1_2 == 1
% THERE IS SOMETHING WRONG IN EITHER THE DATA OR MY LIKELIHOOD. CHECK WITH
% HLW CODE.
settings = struct;
settings.smoother = 1;
settings.covid_dummies  = covid_dummies(startdate:enddate,:);
settings.detrend_y = 1;
settings.states_smoother = 1:3;

% compute intial values of the state and prediction error.
% use HP filter 
xi_0 = zeros(3,1);
xi_0(1:3) = nat_output(startdate-1:-1:startdate-3);

settings.xi_0 = xi_0;
P_0 = 0.2*eye(3);
settings.P_0 = P_0;
settings.fixed_param = [];

% define function to minimize
f_step1 = @(theta) -likelihood_kalman_covid(y_mat,x_mat_step1,theta,@map_theta_to_mats_HLM_step1_covid,settings);

% take starting values from laubach williams
% a_y1 = theta(1);
% a_y2 = theta(2);
% g = theta(3);
% b_y = theta(4);
% b_pi = theta(5);
% sigma_ys = theta(6);
% sigma_yt = theta(7);
% sigma_pi = theta(8);
% phi = theta(9);
% kappa_vec = [1;theta(10:12)];

x0_step1 = [1.5077186; -0.5820255;  0.6793814;  0.1099051;  0.7401697;  0.4137552;  0.7894324;  0.5938689; -0.1093455;  8.5074377;  1.5846599;  1.3925352];
LB_step1 = [-2; -2; 0; 0.025; 0; 0.01; 0.01; 0.01;-inf;1;1;1];
UB_step1 = [2; 2; inf; inf; 1; inf; inf; inf;inf;inf;inf;inf];
[theta_opt_aux,fval,exitflag,~] = fminsearchbnd(f_step1,x0_step1,LB_step1,UB_step1,options);
estimates_step1_aux = evaluate_kalman_covid(y_mat,x_mat_step1,theta_opt_aux,@map_theta_to_mats_HLM_step1_covid,settings);

% Do same procedure as in HLW: they first estimate the full thing, get the
% filtered uncertainty, and rerun everything using that as initial value

settings.P_0 = squeeze(estimates_step1_aux.state_one_ahead_pred(1,:,:));
[theta_opt,fval,exitflag,~] = fminsearchbnd(f_step1,theta_opt_aux,LB_step1,UB_step1,options);

% this is not EXACTLY the same because we are evaluating at slightly
% different parameters. check: use the same parameters, get exactly the same!
%theta_opt_use = [1.4912836 -0.5651727  0.6812395  0.1068789  0.7398945  0.4230496  0.7901214  0.5898836 -0.1093963  8.2716439  1.5898603  1.4021586];
estimates_step1 = evaluate_kalman_covid(y_mat,x_mat_step1,theta_opt,@map_theta_to_mats_HLM_step1_covid,settings);
%% Step 2

% Use estimates in step 1 to compute lambda_g
pot_smoothed_step_1 = estimates_step1.state_smoothed(:,1)/100;
lambda_g = median_unbiased_estimator_stage1(pot_smoothed_step_1);

settings_step2 = struct;
settings_step2.smoother = 1;
settings_step2.detrend_y = 0;
settings_step2.covid_dummies  = covid_dummies(startdate:enddate,:);
settings_step2.states_smoother = 1:6;

xi_0(1:3) = nat_output(startdate-1:-1:startdate-3);
xi_0(4) = nat_output(startdate-1)-nat_output(startdate-2);
xi_0(5) =  nat_output(startdate-2)-nat_output(startdate-3);
xi_0(6) =  nat_output(startdate-3)-nat_output(startdate-4);
settings_step2.xi_0 = xi_0;
settings_step2.P_0 = 0.2*eye(6);
settings_step2.fixed_param = [lambda_g];

% define function to minimize
f_step2 = @(theta) -likelihood_kalman_covid(y_mat,x_mat_step2,theta,@map_theta_to_mats_HLM_step2_covid,settings_step2);

% take starting values from laubach williams
% a_y1 = theta(1);
% a_y2 = theta(2);
% a_r = theta(3);
% a_0 = theta(4);
% a_g = theta(5);
% b_pi = theta(6);
% b_y = theta(7);
% sigma_yt = theta(8);
% sigma_pi = theta(9);
% sigma_ys = theta(10);
% phi = theta(11);
% kappa_vec = [1;theta(12:14)];

x0_step2 = [1.36423315; -0.43058745; -0.08128108; -0.55318179;  0.96902147;  0.68882976;  0.07500347;  0.44800952;  0.79100686;  0.48482401; -0.11274871;
7.35161448;  1.63069704;  2.00573111];
LB_step2 = [-2; -2; -inf; -inf; -inf; 0; .0025; 0; 0; 0;-inf;1;1;1];
UB_step2 = [2; 2; -0.0025; inf; inf; 1; inf; inf; inf; inf;inf;inf;inf;inf];

[theta_opt2_aux,~,~,~] = fminsearchbnd(f_step2,x0_step2,LB_step2,UB_step2,options);
estimates_step2_aux = evaluate_kalman_covid(y_mat,x_mat_step2,theta_opt2_aux,@map_theta_to_mats_HLM_step2_covid,settings_step2);

settings_step2.P_0 = squeeze(estimates_step2_aux.state_one_ahead_pred(1,:,:));
[theta_opt2,fval,exitflag,~] = fminsearchbnd(f_step2,theta_opt2_aux,LB_step2,UB_step2,options);
estimates_step2 = evaluate_kalman_covid(y_mat,x_mat_step2,theta_opt2,@map_theta_to_mats_HLM_step2_covid,settings_step2);

%% Step 3
% Use estimates in step 1 to compute lambda_g

% this is taken from HLW code.
% Y = smoothed output gap.
% X = all regressors in the IS
output_gap_aux = gdp_tot(startdate:enddate)-estimates_step2.state_smoothed(:,1)-theta_opt2(11)*covid_index(startdate:enddate);
y_lambdaz = [output_gap_aux(3:end)];
x_lambdaz = [output_gap_aux(2:end-1),output_gap_aux(1:end-2),...
    (x_mat_step2(3:end,3)+x_mat_step2(3:end,4))/2,...
    4*estimates_step2.state_smoothed(3:end,4), ones(length(y_lambdaz),1)];

% reasonable lambda_z
kappa_vec = settings_step2.covid_dummies(3:end,:)*[1;theta_opt2(end-2:end)]; 
lambda_z = median_unbiased_estimator_stage2_covid(y_lambdaz, x_lambdaz,kappa_vec);

else
    lambda_g = 0.0678583877361746;
    lambda_z = 0.0176786502499355;
end
%%
settings_step3 = struct;
settings_step3.detrend_y = 0;
settings_step3.smoother = 1;
xi_0 = zeros(9,1);
xi_0(1:3) = nat_output(startdate-1:-1:startdate-3);
xi_0(4) = nat_output(startdate-1)-nat_output(startdate-2);
xi_0(5) =  nat_output(startdate-2)-nat_output(startdate-3);
xi_0(6) =  nat_output(startdate-3)-nat_output(startdate-4);
settings_step3.xi_0 = xi_0;
settings_step3.P_0 = 0.2*eye(9);
settings_step3.states_smoother = 1:9;

settings_step3.fixed_param = [lambda_g,lambda_z];
settings_step3.covid_dummies = covid_dummies(startdate:enddate,:);
% define function to minimize
f_step3 = @(theta) -likelihood_kalman_covid(y_mat,x_mat_step3,theta,@map_theta_to_mats_HLM_step3_covid,settings_step3);

% take starting values from laubach williams
% a_y1 = theta(1);
% a_y2 = theta(2);
% a_r = theta(3);
% b_pi = theta(4);
% b_y = theta(5);
% sigma_yt = theta(6);
% sigma_pi = theta(7);
% sigma_ys = theta(8);
% phi = theta(9);
% c = theta(10);
% kappa_vec = [1;theta(11:13)];

%LB_step3 = [-2; -2; -inf; 0; 0.025; 0; 0; 0;-inf;1;1;1;1];
%UB_step3 = [2; 2; -0.0025; 1; inf; inf; inf; inf; inf; inf; inf; inf; inf];

x0_step3 = [1.41452506; -0.48147306; -0.06845139;... % a_y1, a_y2, a_r
      0.68908558;  0.07959243;... %b_pi b_y
      0.43876088;  0.79135821;  0.50258058;... % sigma_yt, sigma_pi, sigma_ys
      -0.11185722;  1.09630223; 7.68332516;  1.58709731 ; 1.96449242]; % phi c kappa_vec

fix_params_step3 = [0; 0; 0;... % a_y1, a_y2, a_r
    0; 0;... %b_pi b_y
    0; 0; 0;... % sigma_yt, sigma_pi, sigma_ys
    0;0;0;0;0]; % phi c kappa_vec


LB_step3 = [-2; -2; -inf;... % a_y1, a_y2, a_r
    0; 0.025;... %b_pi b_y
    0; 0; 0;... % sigma_yt, sigma_pi, sigma_ys
    -inf;0;1;1;1]; % phi c kappa_vec

UB_step3 = [2; 2; -0.0025;... % a_y1, a_y2, a_r
    1; inf;... %b_pi b_y
    inf; inf; inf;... % sigma_yt, sigma_pi, sigma_ys
    inf; inf; inf; inf; inf]; % phi c kappa_vec

fixed_value_step3 = [1.41452506; -0.48147306; -0.06845139;... % a_y1, a_y2, a_r
      0.68908558;  0.07959243;... %b_pi b_y
      0.43876088;  0.79135821;  0.50258058;... % sigma_yt, sigma_pi, sigma_ys
      -0.11185722;  1.09630223; 7.68332516;  1.58709731 ; 1.96449242]; % phi c kappa_vec

% fix parameters by imposing the constraint.
for i_p = 1:length(x0_step3)
    if fix_params_step3(i_p) == 1
    LB_step3(i_p) = fixed_value_step3(i_p);
    UB_step3(i_p) = fixed_value_step3(i_p);
    end
end
 
[theta_opt3_aux,~,~,~] = fminsearchbnd(f_step3,x0_step3,LB_step3,UB_step3,options);
estimates_step3_aux = evaluate_kalman_covid(y_mat,x_mat_step3,theta_opt3_aux,@map_theta_to_mats_HLM_step3_covid,settings_step3);


settings_step3.P_0 = squeeze(estimates_step3_aux.state_one_ahead_pred(1,:,:));
[theta_opt3,fval,exitflag,~] = fminsearchbnd(f_step3,theta_opt3_aux,LB_step3,UB_step3,options);


estimates_step3 = evaluate_kalman_covid(y_mat,x_mat_step3,theta_opt3,@map_theta_to_mats_HLM_step3_covid,settings_step3);

%% Compute parameter SEs and states SEs
N_draws = 1000;
fun_se = @(theta) evaluate_kalman_covid(y_mat,x_mat_step3,theta,@map_theta_to_mats_HLM_step3_covid,settings_step3);
constraints.UB = UB_step3;
constraints.LB = LB_step3;
%a function that is positive if the constraint is NOT satisfied. This is on
%top of the upper and lower bounds.
constraints.const_f = @(theta) (theta(1)+theta(2)>1);
estimates_se = compute_se_kalman(theta_opt3,constraints,fun_se,N_draws);

%% Report Results
clc
fprintf('a_{y,1} = %5.3f \n',theta_opt3(1))
fprintf('a_{y,2} = %5.3f \n',theta_opt3(2))
fprintf('a_{r} = %5.3f \n',theta_opt3(3))
fprintf('b_{pi} = %5.3f \n',theta_opt3(4))
fprintf('b_{y} = %5.3f \n',theta_opt3(5))
fprintf('sigma_{yt} = %5.3f \n',theta_opt3(6)) 
fprintf('sigma_{pi} = %5.3f \n',theta_opt3(7))
fprintf('sigma_{ys} = %5.3f \n',theta_opt3(8)) 
fprintf('phi = %5.3f \n',theta_opt3(9)) 
fprintf('c = %5.3f \n',theta_opt3(10)) 
fprintf('kappa_{2020} = %5.3f \n',theta_opt3(11)) 
fprintf('kappa_{2021} = %5.3f \n',theta_opt3(12)) 
fprintf('kappa_{2022} = %5.3f \n',theta_opt3(13)) 
fprintf('log_likehood = %5.3f \n',fval) 


%% Plot Results

% levels
r_star = 4*theta_opt3(10)*estimates_step3.state_contemporaneous(:,4)+estimates_step3.state_contemporaneous(:,7);
y_natural = estimates_step3.state_contemporaneous(:,1)+theta_opt3(9)*x_mat_step3(:,7); %excluding covid
y_gap = gdp_tot(startdate:enddate)-y_natural;
r_gap = (vardata_aux(startdate:enddate,3)-r_star);

g_trend =  4*theta_opt3(10)*estimates_step3.state_contemporaneous(:,4);

%SE
r_star_se = sqrt(squeeze(estimates_step3.state_contemporaneous_pred(:,4,4))+...
    2*squeeze(estimates_step3.state_contemporaneous_pred(:,4,7))+...
    squeeze(estimates_step3.state_contemporaneous_pred(:,7,7)));

y_gap_se = sqrt(squeeze(estimates_step3.state_contemporaneous_pred(:,1,1)));

dates_use = date(startdate:enddate);
%% Replication of HLW Figure 1: very good!
cd([path vintage task '/_results'])
close all

figure('Position',[0 100 1000 600])

subplot(2,1,1)
plot(dates_use,y_gap,'LineWidth',3,'Color','k')
yline(0,LineStyle="--")
hold on
if recent_data ~=3
plot(dates_use,r_gap,'LineWidth',3,'Color','b','LineStyle','--')
rr_str = 'Real rate Gap';
else
    rr_str = '';
end
legend({'Output Gap',rr_str},'Location','northwest')
title('Output Gap and Real Rate Gap','Interpreter','latex')

subplot(2,1,2)
plot(dates_use,r_star,'LineWidth',3,'Color','k')
hold on
plot(dates_use,g_trend ,'LineWidth',3,'Color','b','LineStyle','--')
yline(0,LineStyle="--")
legend({'$r^{*}$','Trend Growth'},'Interpreter','latex','Location','southwest')
title('$r^{*}$ and Trend Growth','Interpreter','latex')

if save_fig == 1
    print(['replication_HLW_' num2str(floor(date(startdate))) '_' num2str(floor(date(enddate)))],'-dpng')
    print(['replication_HLW_' num2str(floor(date(startdate))) '_' num2str(floor(date(enddate)))],'-depsc')
end

% figure(1)
% plot(dates_use,y_gap,'LineWidth',3,'Color','k')
% hold on
% plot(dates_use,r_gap,'LineWidth',3,'Color','b','LineStyle','--')
% yline(0,LineStyle="--")
% legend({'Output Gap','Real rate Gap'},'Location','northwest')
% title('Output Gap and Real Rate Gap','Interpreter','latex')
% 
% if save_fig == 1
%     print('replication_HLW_gaps','-dpng')
%     print('replication_HLW_gaps','-depsc')
% end
% 
% figure(2)
% plot(dates_use,r_star,'LineWidth',3,'Color','k')
% hold on
% plot(dates_use,g_trend ,'LineWidth',3,'Color','b','LineStyle','--')
% yline(0,LineStyle="--")
% legend({'$r^{*}$','Trend Growth'},'Interpreter','latex')
% title('$r^{*}$ and Trend Growth','Interpreter','latex')
% 
% if save_fig == 1
%     print('replication_HLW_r_star','-dpng')
%     print('replication_HLW_r_star','-depsc')
% end
