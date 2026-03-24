%This code constructs fci* and r* directly from y_gap using the dynamic output-asset price equations
% IS in HLW, OA in CCS. 

% Result: errors are playing a big role in both cases.


%% HOUSEKEEPING
 
clc, clear all, close all

warning('off','MATLAB:dispatcher:nameConflict')
path = '/Users/tomyc/Dropbox (MIT)/FCI_targeting/FCIstar/code/empirics';
vintage = '/25_06_14';
task = '\construct_fcistar_ygap';

data_path = '/Users/tomyc/Dropbox (MIT)/FCI_targeting/code/empirics';

addpath([path vintage '/_aux_fun'])
addpath([path '/_data'])
addpath([path '/_data/HLW_data'])
cd([path vintage task])
%% Settings
save_fig = 0;
recent_data = 1;
run_step_1 = 0; % 1 if you want to rerun step 1, 0 if you are fine calibrating the lambda_g
backward_inf = 1; %use the 1 year past inflation as proxt for future inflation.

set_covid_2022 = 1; %1 sets the covid dummy to 0 after 2022.

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

% Load HLW
T_HLW_y_gap = readtable('Holston_Laubach_Williams_current_estimates.xlsx','Sheet','HLW Estimates');

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


HLW_dates = T_HLW_y_gap{:,1};
startdate_HLW = find(HLW_dates==startdate_HLW_str);
enddate_HLW = find(HLW_dates==enddate_HLW_str);
HLW_y_gap = T_HLW_y_gap{startdate_HLW:enddate_HLW,15};
HLW_r_star_aux =  T_HLW_y_gap{startdate_HLW-2:enddate_HLW,11};
HLW_r_star_smooth = (HLW_r_star_aux(2:end-1)+HLW_r_star_aux(1:end-2))/2;

dates_use = date(startdate:enddate);
[nat_output,output_gap_hp] =hpfilter(gdp_tot,36000); % separate in "natural" and "gap"
%start p_t such that p_t - p^* is close to zero.
%p_0 = nat_output(startdate-4+n_lags)+(nat_output(startdate-4+n_lags)-nat_output(startdate-5+n_lags))-2*log(0.02);
fci_use = fci_2;
vardata_aux = [y_gap_cbo fci_use];
vardata = vardata_aux(startdate:enddate,:);
vardata_2_aux = [y_gap_cbo r_rate];
vardata_2 = vardata_2_aux(startdate-2:enddate,:);

a_eta = 0.26;


% construct fci star
fci_star =vardata(2:end,2) + (vardata(2:end,1)-vardata(1:end-1,1))/a_eta;

% construct r_star using the dynanmic IS from HLW
a_1 = 1.417;
a_2 = -0.483;
a_r = -0.068; 

r_star_avg = (vardata_2(2:end-1,2)+vardata_2(1:end-2,2))/2 - ...
    (vardata_2(3:end,1) - a_1 * vardata_2(2:end-1,1) - a_2* vardata_2(1:end-2,1))/a_r; 

%%
% plot

figure
plot(dates_use(2:end),fci_star);
hold on
plot(dates_use(2:end),vardata(2:end,2));
%%
figure
plot(dates_use,HLW_r_star_smooth );
hold on
plot(dates_use,r_star_avg);

% main result here is that this does not work well for either us or HLW.
