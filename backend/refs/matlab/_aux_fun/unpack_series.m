function  out = unpack_series(estimates_step3,x_mat_step3,covid_d_pos_step3,fci_pos,y_n_pos,g_pos,delta_pos,startdate,enddate,gdp_tot,param)

fci_star_true =  estimates_step3.state_contemporaneous(:,fci_pos);
% measured FCI star
fci_star = fci_star_true;
y_star= estimates_step3.state_contemporaneous(:,y_n_pos+1)+estimates_step3.state_contemporaneous(:,g_pos+1)...
    + estimates_step3.state_contemporaneous(:,delta_pos)-param.rho_delta*estimates_step3.state_contemporaneous(:,delta_pos+1)...
    +param.phi*x_mat_step3(:,covid_d_pos_step3(1)); %including covid
y_gap = gdp_tot(startdate:enddate)-y_star;
% y_star_smoothed = estimates_step3.state_smoothed(:,y_n_pos+1)+estimates_step3.state_smoothed(:,g_pos+1)...
%     + estimates_step3.state_smoothed(:,delta_pos)-param.rho_delta*estimates_step3.state_smoothed(:,delta_pos+1)...
%     +param.phi*x_mat_step3(:,covid_d_pos_step3(1)); %including covid
% y_gap_smoothed = gdp_tot(startdate:enddate) -y_star_smoothed;
% % y_gap_smoothed = gdp_tot(startdate:enddate)-y_star_smoothed;
% fci_gap = fci_use(startdate:enddate)-fci_star;

out = struct;

out.y_gap = y_gap;
out.fci_star = fci_star;

end