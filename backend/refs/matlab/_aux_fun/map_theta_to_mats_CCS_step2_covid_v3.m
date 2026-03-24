function matrices = map_theta_to_mats_CCS_step2_covid_v3(theta,fixed_param,covid_dummies)

% Write step 1 state space system following replication code in HLM.

% unpack parameters
a_y1 = theta(1);
a_y2 = theta(2);
a_p = theta(3);
a_0 = theta(4);
a_g = theta(5);
b_pi = theta(6);
b_y = theta(7);
sigma_yt = theta(8);
sigma_pi = theta(9);
sigma_ys = theta(10);
phi = theta(11);
kappa_vec = [1;theta(12:14)];
sigma_pi_t = theta(15);

lambda_g = fixed_param(1);

F = [1 0 0 1 0 0;
     1 0 0 0 0 0;
     0 1 0 0 0 0;
     0 0 0 1 0 0;
     0 0 0 1 0 0;
     0 0 0 0 0 1];

Q = [sigma_ys^2 0 0 0 0 0;
            0 0 0 0 0 0;
            0 0 0 0 0 0;
            0 0 0 (lambda_g*sigma_ys)^2 0 0;
            0 0 0 0 0 0;
            0 0 0 0 0 sigma_pi_t^2];

H_t = [1    -a_y1     -a_y2 0 a_g 0;
        0       -b_y        0     0   0 (1-b_pi)];

H = H_t';

A_t = [a_y1 a_y2 a_p    0     a_0 phi -phi*a_y1 -phi*a_y2;
       b_y   0      0   b_pi  0 0 -phi*b_y 0];

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