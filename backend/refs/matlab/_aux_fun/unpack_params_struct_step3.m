function param = unpack_params_struct_step3(theta)

param = struct;

% unpack parameters
param.a_y1 = theta(1);
param.a_y2 = theta(2);

param.b_pi = theta(3);
param.b_y = theta(4);
param.sigma_yt = theta(5);
param.sigma_pi = theta(6);
param.sigma_ys = theta(7);
param.sigma_pi_e = theta(8);
param.sigma_delta = theta(9);
param.phi = theta(10);
param.kappa_vec = [1;theta(11:13)];
param.rho_pi = theta(14);
param.rho_delta = theta(15);
param.eta = theta(16);
param.a_f = -theta(17)*(1-theta(16))/(1-theta(16)^4);
param.alpha = 1;
param.c_g = 1;
end