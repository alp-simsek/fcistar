function matrices = map_theta_to_mats_CCS_struct_step1(theta,fixed_param,covid_dummies)

% Write step 1 state space system following replication code in HLM.

param =  unpack_params_struct_step1(theta);

% unpack parameters
a_y1 = param.a_y1;
a_y2 = param.a_y2;
b_pi = param.b_pi;
b_y = param.b_y;
g = param.g;
sigma_yt = param.sigma_yt; % signal
sigma_ys = param.sigma_ys;
sigma_pi = param.sigma_pi;
sigma_pi_e = param.sigma_pi_e;
phi = param.phi;
kappa_vec = param.kappa_vec;
rho_pi = param.rho_pi;
c_pi = 2;

F = zeros(6,6);
F(1,1) = 1; F(2,1)=1; F(3,2) = 1; F(4,3) = 1;
F(5,5) = rho_pi;F(5,6)=c_pi*(1-rho_pi); F(6,6) = 1;


Q = zeros(6,6);
Q(1,1) = sigma_ys^2;
Q(5,5) = sigma_pi_e^2;


H_t = [1 -(1+a_y1) -(a_y2-a_y1) a_y2    0      0;
        0 -b_y          0     0    (1-b_pi)    0 ];

H = H_t';

A_t = [1+a_y1 a_y2-a_y1 -a_y2     0    phi -phi*(1+a_y1) -phi*(a_y2-a_y1) phi*a_y2;
        b_y       0        0    b_pi    0   -phi*b_y          0               0];

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