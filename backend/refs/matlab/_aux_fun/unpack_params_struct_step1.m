function param = unpack_params_struct_step1(theta)

param = struct;

% unpack parameters
param.a_y1 = theta(1);
param.a_y2 = theta(2);
param.b_pi = theta(3);
param.b_y = theta(4);
param.g = theta(5);
param.sigma_yt = theta(6);
param.sigma_pi = theta(7);
param.sigma_ys = theta(8);
param.sigma_pi_e = theta(9);
param.phi = theta(10);
param.kappa_vec = [1;theta(11:13)];
param.rho_pi = theta(14);

end