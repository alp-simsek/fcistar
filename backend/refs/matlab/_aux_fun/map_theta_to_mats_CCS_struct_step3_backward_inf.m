function matrices = map_theta_to_mats_CCS_struct_step3_backward_inf(theta,fixed_param,covid_dummies)

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
c_g = 1;

lambda_g = fixed_param(1);
sigma_g = lambda_g*sigma_ys;


% states order:
% xi_t = [y_{t}^{*},y_{t-1}^{*},y_{t-2}^{*},g_{t},g_{t-1},g_{t-2},FCI_{t}^{*},FCI_{t-1}^{*},\delta_{t},\delta_{t-1},\delta_{t-2}]


F = zeros(11,11);
F(1,1) = 1; F(1,4) = 1; % y* equation
F(2,1)=1; F(3,2) = 1; % y* keep track of the lags
F(4,4) = 1; % g equation
F(5,4) = 1; % keep track of g lags
F(6,5) = 1;

F(9,9) = rho_delta; % delta eq
F(10,9) = 1; %keep track of delta lags
F(11,10) = 1;

% coefficeints on FCI* equation. In the code a_f is defined with a minus,
% so change the signs with respect to the paper.

F(7,1) = -eta/a_f; F(7,2) = eta/a_f; % coefficients on y*
F(7,4) = c_g*(1-eta)/a_f; F(7,5) = c_g*eta/a_f; % coefficients on g
F(7,7) = eta;
F(7,9) = (1/a_f)*(eta+(1-rho_delta)*rho_delta); F(7,10) = -(1/a_f)*eta*rho_delta;
F(8,7) = 1; % FCI* keep track of the lags


Q = zeros(11,11);
Q(1,1) = sigma_ys^2; % y* equation
Q(4,4) = sigma_g^2; % g equation
Q(9,9) = sigma_delta^2;  % delta eq

%covariances of  FCI* innovations
Q(7,1) = (1/ a_f)*sigma_ys^2;
Q(1,7) = Q(7,1);
Q(7,4) = c_g*(1/ a_f)*sigma_g^2;
Q(4,7) = Q(7,4);
Q(7,9) = -(1/ a_f)*(eta+rho_delta)*sigma_delta^2;
Q(9,7) = Q(7,9);
Q(7,7) = (1/ a_f)^2*((c_g)^2*sigma_g^2+sigma_ys^2+((eta+rho_delta)^2)*sigma_delta^2); % FCI star eq


H_t = zeros(2,11);
H_t(1,2) = 1; H_t(1,3) = -1; %coefficeints on natural output
H_t(1,5) = 1; H_t(1,6) = -1; %coefficients on lagged growth
H_t(1,7) = - a_f; %coefficients on FCI*
H_t(1,9) = 1; %coefficient on delta_t
H_t(1,10) = -(1+rho_delta); %coefficient on delta_t-1
H_t(1,11) = rho_delta; %coefficient on delta_t-2

% H_t(2,2) = -b_y; %make only natural outupt show up in the PC
H_t(2,3) = -b_y; %Phillips curve
H_t(2,6) = -b_y;
H_t(2,10) = -b_y;
H_t(2,11) = b_y*rho_delta;

H = H_t';

% x_t = [y_{t-1},FCI_{t-1}^{m},\pi_{t-1},\pi_{t-2:4},d_{t},d_{t-1}]

A_t = zeros(2,6);
A_t(1,1) = 1; A_t(1,2) = a_f;
A_t(1,5) = phi; A_t(1,6) = -phi;
A_t(2,1) = b_y; A_t(2,3) = b_pi; A_t(2,4) = (1-b_pi); A_t(2,6) = -b_y*phi;

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