function matrices = map_theta_to_mats_CCS_struct_step3(theta,fixed_param,covid_dummies)

% Write step 1 state space system following replication code in HLM.

param =  unpack_params_struct_step3(theta);

% unpack parameters
a_y1 = param.a_y1;
a_y2 = param.a_y2;
a_y3 = a_y2;
a_f = param.a_f;
b_pi = param.b_pi;
b_y = param.b_y;
sigma_yt = param.sigma_yt; % signal
sigma_ys = param.sigma_ys;
sigma_pi = param.sigma_pi;
sigma_pi_e = param.sigma_pi_e;
sigma_delta = param.sigma_delta;
phi = param.phi;
kappa_vec = param.kappa_vec;
rho_delta = param.rho_delta;
rho_pi = param.rho_pi;
eta = param.eta;
c_pi = param.c_pi;


alpha = 1;
lambda_g = fixed_param(1);
sigma_g = lambda_g*sigma_ys;

% states order:
% [y_{t}^{*},y_{t-1}^{*},y_{t-2}^{*},y_{t-3}^{*},g_{t},FCI_{t}^{*},FCI_{t-1}^{*},FCI_{t-2}^{*},FCI_{t-3}^{*},\delta_{t},\pi_{t}^{e},1]

F = zeros(12,12);
F(1,1) = 1; F(1,5) = 1; % y* equation
F(2,1)=1; F(3,2) = 1; F(4,3) = 1; % y* keep track of the lags
F(5,5) = 1; % g equation
F(6,5) = (1-eta)/a_f; F(6,6) = eta; F(6,10) = (1-rho_delta)*rho_delta/a_f; % FCI* eq
F(7,6) = 1; F(8,7) = 1; F(9,8) = 1; % FCI* keep track of the lags
F(10,10) = rho_delta; % delta eq
F(11,11) = rho_pi; F(11,12) =(1-rho_pi)*c_pi;   % pi_e eq 
F(12,12) = 1;




Q = zeros(12,12);
Q(1,1) = sigma_ys^2; % y* equation
Q(5,5) = sigma_g^2; % g equation
Q(10,10) = sigma_delta^2;  % delta eq
Q(11,11) = sigma_pi_e^2; %sigma pi_e eq

Q(6,6) = (1/ a_f)^2*(sigma_g^2+eta^2*sigma_ys^2+(1-rho_delta)^2*sigma_delta^2); % FCI star eq
Q(6,1) = -(1/ a_f)*eta*sigma_ys^2; % cov with y*
Q(1,6) = Q(6,1);
Q(6,5) = (1/ a_f)*sigma_g^2; % cov with g
Q(5,6) = Q(6,5);
Q(6,10) = (1/ a_f)*(1-rho_delta)^2*sigma_delta^2;  % cov with delta
Q(10,6) = Q(6,10);

H_t = zeros(2,12);
H_t(1,1) = 1; H_t(1,2) = -(1+a_y1); H_t(1,3) = -(a_y2-a_y1); H_t(1,4) = a_y3; %coefficeints on natural output
H_t(1,7) = - alpha*a_f; H_t(1,8) = alpha*a_f*a_y1; H_t(1,9) =  alpha*a_f*a_y3; %coefficients on FCI*
H_t(2,2) = -b_y; H_t(2,11) = (1-b_pi); %Phillips curve

H = H_t';

% x_{t}=[y_{t-1},y_{t-2},y_{t-3},FCI_{t-1}^{m},\pi_{t-1},d_{t},d_{t-1},d_{t-2},d_{t-3}]

A_t = [1+a_y1 (a_y2-a_y1) -a_y3     alpha*a_f 0     phi -phi*(1+a_y1) -phi*(a_y2-a_y1) phi*a_y3;
        b_y       0        0            0        b_pi    0   -phi*b_y          0               0];

A = A_t';


R =  [sigma_yt^2 0;
            0   sigma_pi^2];

matrices.F = F;
matrices.A = A;
matrices.H = H;
matrices.Q = Q;
matrices.R = R;
matrices.kappa_seq = covid_dummies*kappa_vec; % this is the vector of kappas to factor the R_t

end