function matrices = map_theta_to_mats_CCS_step3_covid_v2(theta,fixed_param,covid_dummies)

% Write step 1 state space system following replication code in HLM.

% unpack parameters
a_y1 = theta(1);
a_y2 = theta(2);
a_p = theta(3);
b_pi = theta(4);
b_y = theta(5);
sigma_yt = theta(6);
sigma_pi = theta(7);
sigma_ys = theta(8);
phi = theta(9);
c_p = theta(10);
kappa_vec = [1;theta(11:13)];
sigma_pi_t = theta(14);
rho_ps = 1;


lambda_g = fixed_param(1);
lambda_z = fixed_param(2);

F = [1 0 0 1 0 0 0 0 0;
     1 0 0 0 0 0 0 0 0;
     0 1 0 0 0 0 0 0 0;
     0 0 0 1 0 0 0 0 0;
     0 0 0 1 0 0 0 0 0;
     0 0 0 0 0 rho_ps 0 0 0;
     0 0 0 0 0 1 0 0 0;
     0 0 0 0 0 0 0 1 0;
     0 0 0 0 0 0 0 1 0;];

Q = zeros(9,9);
Q(1,1) = sigma_ys^2;
Q(4,4) = (lambda_g*sigma_ys)^2;
Q(6,6) = (lambda_z*sigma_yt/(c_p*a_p))^2;
Q(8,8) = sigma_pi_t^2;

A_t = [a_y1 a_y2    a_p     0       phi -phi*a_y1 -phi*a_y2;
       b_y   0      0     b_pi    0 -phi*b_y 0];

A = A_t';

H_t = [1  -a_y1 -a_y2   0 -a_p*c_p 0 -a_p   0    0;
        0   -b_y    0   0      0   0        0  1 -b_pi];

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