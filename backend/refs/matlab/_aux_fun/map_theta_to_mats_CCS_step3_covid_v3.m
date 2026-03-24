function matrices = map_theta_to_mats_CCS_step3_covid_v3(theta,fixed_param,covid_dummies)

% Write step 1 state space system following replication code in HLM.

lambda_g = fixed_param(1);

params = unpack_params(theta);

% unpack parameters
a_y1 = params.a_y1;
a_y2 = params.a_y2;
a_p = params.a_p;
b_pi = params.b_pi;
b_y = params.b_y;
sigma_yt = params.sigma_yt;
sigma_pi = params.sigma_pi;
sigma_ys = params.sigma_ys;
sigma_pi_t = params.sigma_pi_t;
sigma_x = params.sigma_x;
phi = params.phi;
kappa_vec = params.kappa_vec;
rho_f1 = params.rho_f1;
rho_f2 = params.rho_f2;
c_xg =  params.c_xg; %coefficient of eps_g on f^* equation
c_xys =  params.c_xys ; %coefficient of eps_y* on f^* equation
c_e_pi =  params.c_e_pi;
rho_pi =  params.rho_pi;
c_pi = 2;


% states:
% [y*,y*_{-1},y*_{-2},g,g_{-1},f^*,f^*_{-1},pi^{e},1]


F = zeros(9,9);
F(1,1) = 1; F(1,4) = 1;
F(2,1) = 1;
F(3,2) = 1;
F(4,4) = 1;
F(5,4) = 1;
F(6,6) = rho_f1;
F(6,7) = rho_f2;
F(7,6) = 1;
F(8,8) = rho_pi; F(8,9) = c_pi*(1-rho_pi);
F(9,9) = 1;
% 
% F = [1 0 0 1 0 0     0   0  0;
%      1 0 0 0 0 0     0   0  0;
%      0 1 0 0 0 0     0   0  0;
%      0 0 0 1 0 0     0   0  0;
%      0 0 0 1 0 0     0   0  0;
%      0 0 0 0 0 rho_x1 rho_x2 0 0;
%      0 0 0 0 0 1     0   0  0;
%      0 0 0 0 0 0     0   rho_pi c_pi*(1-rho_pi);
%      0 0 0 0 0 0     0   0 0 1];



Q = zeros(9,9);
Q(1,1) = sigma_ys^2;
Q(4,4) = (lambda_g*sigma_ys)^2;

% variance in the f* equation
Q(6,6) = sigma_x^2 + (c_xg*lambda_g*sigma_ys)^2+ (c_xys*sigma_ys)^2 + (c_e_pi*sigma_pi_t)^2;

% cov with growth shocks
Q(4,6) = c_xg*(lambda_g*sigma_ys)^2;
Q(6,4) = Q(4,6);

% cov with natural output (supply) shocks
Q(1,6) = c_xys*(sigma_ys^2);
Q(6,1) = Q(1,6);

%cov with expected inflation shocks
Q(6,8) = c_e_pi*(sigma_pi_t)^2;
Q(8,6) = c_e_pi*(sigma_pi_t)^2;
Q(8,8) = sigma_pi_t^2;



A_t = [a_y1 a_y2    a_p     0       phi -phi*a_y1 -phi*a_y2;
       b_y   0      0     b_pi    0 -phi*b_y 0];

A = A_t';


H_t = [1  -a_y1 -a_y2   0 0 0      -a_p   0      0;
        0   -b_y  0   0      0   0        0  (1-b_pi) 0];

H = H_t';

R =  [sigma_yt^2 0;
            0   sigma_pi^2];

matrices.F = F;
matrices.A = A;
matrices.H = H;
matrices.Q = Q;
matrices.R = R;
matrices.kappa_seq = covid_dummies*kappa_vec; % this is the vector of kappas to factor the R_t

end