clear all
close all
clc

warning('off','MATLAB:dispatcher:nameConflict')

path = '/Users/tomyc/Dropbox (MIT)/FCI_targeting/FCIstar/code/empirics';
vintage = '/26_01_10';
task = '/assemble_data';
data_path_save = [path '/_data'];
%data_path_fci = ['/Users/tomyc/Dropbox (MIT)/FCI_targeting/code/empirics/_data/new_fci_index'];
data_path_fci = [data_path_save '/fci_index'];
addpath(data_path_save);
addpath([data_path_save '/covid_stringency_index']);
addpath([data_path_save '/HLW_data']);
addpath(data_path_fci);
cd([path vintage task]);

download_data = 1;

%% Download latest data on inflation, interest rates, output and consumption from Fred.

startdate = '01/01/1959';
enddate = '09/30/2025';
dates_aux = (1959:0.25:2025.5)';

series = {'GDPC1',... % real GDP
            'PCEPILFE',... %core PCE
            'FEDFUNDS',... %Fed Funds rate
            'INTDSRUSM193N',... % FRBNY discount rate 
            'NFCI',...%Chicago Fed financial conditions index 
            'GDPPOT'}; % Potential GDP, CBO

%%

if download_data ==1
url = 'https://fred.stlouisfed.org/';
c = fred(url);
c.DataReturnFormat = 'timetable';
var_labels = cell(length(series),1);
new_data = [];

% get raw data

for i=1:length(series)
d = fetch(c,series{i},startdate, enddate);
% figure out the frequency:
arr_aux = table2array(d.Data{1});
dates_aux_2 = d.Data{1}.Time;

% figure out the frequency
time_interval = days(dates_aux_2(2)-dates_aux_2(1));

if (time_interval<95) & (time_interval>85);
   new_data = [new_data arr_aux];
else %convert to quarterly 
    TT2 = convert2quarterly(d.Data{1},'Aggregation','mean');
    if i==1 || size(TT2,1) == size(new_data,1)
        new_data = [new_data TT2.Var2];
    else
        [year_aux_1,month_aux_1,~] = ymd(TT2.Time(1));
        [year_aux_2,month_aux_2,~] = ymd(TT2.Time(end));

        date_index_1 = year_aux_1+0.25*floor((month_aux_1-1)/3);
        date_index_2 = year_aux_2+0.25*floor((month_aux_2-1)/3);

        position_index_1 = find(dates_aux==date_index_1);
        position_index_2 = find(dates_aux==date_index_2);

        aux_vec = nan(size(new_data,1),1);
        aux_vec(position_index_1:position_index_2) = TT2.Var2;
        new_data = [new_data aux_vec];
    end
end
var_labels{i,1} = d.SeriesID;
end

cd(data_path_save);
save data_fred_raw new_data var_labels
cd([path vintage task]);
else
    load data_fred_raw
end


% if download_data ==1
% url = 'https://fred.stlouisfed.org/';
% c = fred(url);
% c.DataReturnFormat = 'timetable';
% var_labels = cell(length(series),1);
% new_data = [];
% 
% % get raw data
% 
% for i=1:length(series)
% d = fetch(c,series{i},startdate, enddate);
% if strcmp(d.Frequency{1,1},'Quarterly')
%    new_data = [new_data d.Data(:,2)];
% else %convert to quarterly taking end of period
%     TT2 = convert2quarterly(d.Data{1,1},'Aggregation','mean');
%     if i==1 || size(TT2,1) == size(new_data,1)
%         new_data = [new_data TT2.Var1];
%     else
%         [year_aux_1,month_aux_1,~] = ymd(TT2.Time(1));
%         [year_aux_2,month_aux_2,~] = ymd(TT2.Time(end));
% 
%         date_index_1 = year_aux_1+0.25*floor((month_aux_1-1)/3);
%         date_index_2 = year_aux_2+0.25*floor((month_aux_2-1)/3);
% 
%         position_index_1 = find(dates_aux==date_index_1);
%         position_index_2 = find(dates_aux==date_index_2);
% 
%         aux_vec = nan(size(new_data,1),1);
%         aux_vec(position_index_1:position_index_2) = TT2.Var1;
%         new_data = [new_data aux_vec];
%     end
% end
% var_labels{i,1} = d.SeriesID;
% end
% 
% cd(data_path_save);
% save data_fred_raw new_data var_labels
% cd([path vintage task]);
% else
%     load data_fred_raw
% end
%% Load COVID stringency data
% T_covid = readtable('OxCGRT_timeseries_StringencyIndex_v1.csv');
% % row 386 is USA
% 
% t1 = datetime(2020,01,1,1,0,0);
% t2 = datetime(2023,02,28,1,0,0);
% dates_covid = t1:t2; %dates
% usa_daily_covid = T_covid{386,8:end}; %index values
% usa_covid_tt = timetable(dates_covid',usa_daily_covid');
% 
% %convert to quarterly
% usa_covid_tt_q = convert2quarterly(usa_covid_tt,Aggregation=["mean"]);

%% Load COVID data from HLW. After 2024Q4 set it to zero.
TT_HLM = readtable("Holston_Laubach_Williams_current_estimates.xlsx",'Sheet','US input data');
covid_index = zeros(length(dates_aux),1);
startdate_covidhlw = find(dates_aux==1960);
enddate_covid_interp = find(dates_aux==dates_aux(end));
% start from 5th since this data starts in 1960, not 1959.
covid_index(startdate_covidhlw:enddate_covid_interp) = TT_HLM{:,6};

%% Load FED's New FCI 
data_table_new_fci = readtimetable([data_path_fci '/fci_g_public_quarterly_3yr']);
% convert to quarterly

data_fci = table2array(data_table_new_fci);

[year_aux_1,month_aux_1,~] = ymd(data_table_new_fci.date(1));
[year_aux_2,month_aux_2,~] = ymd(data_table_new_fci.date(end));

date_index_fci_1 = year_aux_1+0.25*floor((month_aux_1-1)/3);
date_index_fci_2 = year_aux_2+0.25*floor((month_aux_2-1)/3);

position_index_1 = find(dates_aux==date_index_fci_1);
position_index_2 = find(dates_aux==date_index_fci_2);


dates_fci = (date_index_fci_1:0.25:(date_index_fci_2))';

init_fci = find(dates_aux  == dates_fci(1));
end_fci  =  find(dates_fci == dates_aux(end));

fci_2 = nan(length(dates_aux),1);
fci_2(init_fci:end) = data_fci(1:end_fci);
%%
% also load the 1 year one

data_table_new_fci_1y = readtimetable([data_path_fci '/fci_g_public_quarterly_1yr']);
% convert to quarterly

data_fci_1y = table2array(data_table_new_fci_1y);

[year_aux_1_1y,month_aux_1_1y,~] = ymd(data_table_new_fci_1y.date(1));
[year_aux_2_1y,month_aux_2_1y,~] = ymd(data_table_new_fci_1y.date(end));

date_index_fci_1_1y = year_aux_1_1y+0.25*floor((month_aux_1_1y-1)/3);
date_index_fci_2_1y = year_aux_2_1y+0.25*floor((month_aux_2_1y-1)/3);

position_index_1_1y = find(dates_aux==date_index_fci_1_1y);
position_index_2_1y = find(dates_aux==date_index_fci_2_1y);


dates_fci_1y = (date_index_fci_1_1y:0.25:(date_index_fci_2_1y))';

init_fci_1y = find(dates_aux  == dates_fci_1y(1));
end_fci_1y  =  find(dates_fci_1y == dates_aux(end) );

fci_2_1y = nan(length(dates_aux),1);
fci_2_1y(init_fci:end) = data_fci_1y(1:end_fci);
%%

startdate_covid = find(dates_aux==2020);
enddate_covid = find(dates_aux==2022.75);


gdp_tot = log(new_data(:,1));
pce_index = log(new_data(:,2));
ffr_base = new_data(:,3);
ffr_ny = new_data(:,4);
fci = new_data(:,5);
gdp_pot_cbo = log(new_data(:,6));
% covid_index = zeros(length(gdp_tot),1);
% covid_index(startdate_covid:enddate_covid) = usa_covid_tt_q{1:end-1,1};
% % assume linear decline between 2022.75 and 2024.75;
% linear_decline = covid_index(enddate_covid)*(1:-1/(enddate_covid_interp-enddate_covid):0);
% covid_index(enddate_covid:enddate_covid_interp) = linear_decline';
% express discount rates in 365 days

ffr_base_d = 100*((1+ffr_base/36000).^365 -1);
ffr_ny_d = 100*((1+ffr_ny/36000).^365 -1);
% add also a vector of covid indicator dummies for the variance.



ffr = ffr_base_d;
date_end_ny = find(dates_aux == 1964.75);
ffr(1:date_end_ny) = ffr_ny_d(1:date_end_ny);
%pce_core = [NaN;400*(pce_index(2:end)-pce_index(1:end-1))];
pce_core = [NaN;400*(log(new_data(2:end,2)./new_data(1:end-1,2)))];
expected_pce = movmean(pce_core,[3 0],'includenan');

covid_dummies = zeros(length(gdp_tot),4);
covid_dummies(:,1) = ones(length(gdp_tot),1);
% 2020Q2 to 2020Q4
covid_dummies(startdate_covid+1:startdate_covid+3,:) = repmat([0 1 0 0],3,1);
% 2021Q1 to 2021Q4
covid_dummies(startdate_covid+4:startdate_covid+7,:) = repmat([0 0 1 0],4,1);
% 2022Q1 to 2022Q4
covid_dummies(startdate_covid+8:startdate_covid+11,:) = repmat([0 0 0 1],4,1);


%%
data_2 = [dates_aux gdp_tot pce_core expected_pce ffr covid_index fci fci_2 fci_2_1y gdp_pot_cbo];
covid_dummies_data = [dates_aux covid_dummies];

series_name_2 = {'Date','Real GDP','PCE Core Inflation','Expected Inflation', 'FFR', 'Covid Index',...
    'FCI - Chicago Fed', 'New FCI - 3 year lookback', 'New FCI - 1 years lookback','CBO Potential GDP'};

%% Save VAR Data

data_var = data_2;
series_name_var = series_name_2;

cd(data_path_save)
save data_hlm data_var series_name_var covid_dummies_data







