function [xi_0,P_0] = set_init_val(theta,xi_1,P_1)

% theta is parameter vector, xi_1 is the desired initial state vector for
% [y*,y*_{t-1},y*_{t-2},g_{t},g_{t-1},f^*_{t},pi^{e}_{t}] and P_1 is the
% desired var cov mat.

% main relation is that f^* = theta(10) g + x
% so E_0[x] = E_0[f^*]-theta(10) E_0[g]
% and Var_0[x] = Var0[f^*]

xi_0 = xi_1;
xi_0(6:7) = xi_1(6:7)-theta(10)*xi_1(4:5);

P_0 = P_1;

P_0(6,6) = P_1(6,6)-theta(10)^2*P_1(4,4);
P_0(7,7) = P_1(7,7)-theta(10)^2*P_1(5,5);
P_0(6,7) = theta(16)*P_0(6,6);
P_0(7,6) = theta(16)*P_0(6,6);

P_0(4,6) = theta(10)*P_1(4,4);
P_0(6,4) = P_0(6,4);
end