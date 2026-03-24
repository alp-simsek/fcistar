function param = unpack_params(theta)

param = struct;

% unpack parameters
param.a_y1 = theta(1);
param.a_y2 = theta(2);
param.a_p = theta(3);
param.b_pi = theta(4);
param.b_y = theta(5);
param.sigma_yt = theta(6);
param.sigma_pi = theta(7);
param.sigma_ys = theta(8);
param.sigma_pi_t = theta(9);
param.sigma_x = theta(10);
param.phi = theta(11);
param.kappa_vec = [1;theta(12:14)];
param.rho_f1 = theta(15); % coefficients on AR(2) for f*  
param.rho_f2 = theta(16); % coefficients on AR(2) for f* 
param.c_xg = theta(17); %coefficient of eps_g on f^* equation
param.c_xys = theta(18); %coefficient of eps_y* on f^* equation
param.c_e_pi = theta(19); %coefficient of eps_pi^{e} on f^* equation
param.rho_pi = theta(20); % coefficient on AR(1) ofr pi
param.c_pi = 2;


end