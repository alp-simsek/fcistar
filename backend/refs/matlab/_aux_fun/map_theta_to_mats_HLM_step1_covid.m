function matrices = map_theta_to_mats_HLM_step1_covid(theta,fixed_param,covid_dummies)

% Write step 1 state space system following replication code in HLM.

% unpack parameters
a_y1 = theta(1);
a_y2 = theta(2);
b_pi = theta(3);
b_y = theta(4);
g = theta(5);
sigma_ys = theta(8); % signal
sigma_yt = theta(6);
sigma_pi = theta(7);
phi = theta(9);
kappa_vec = [1;theta(10:12)];

% add one more state, g, and make it constant.
F = [1 0 0;
     1 0 0;
     0 1 0];

Q = [sigma_ys^2 0 0;
            0 0 0 ;
            0 0 0];

H_t = [1 -a_y1 -a_y2;
        0 -b_y    0 ];

H = H_t';

A_t = [a_y1 a_y2 0 0 phi -phi*a_y1 -phi*a_y2;
    b_y 0 b_pi (1-b_pi) 0 -phi*b_y 0];

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