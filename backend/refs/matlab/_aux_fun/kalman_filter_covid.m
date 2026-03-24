function  output = kalman_filter_covid(y,x,F,A,H,Q,R_base,kappa_seq,settings)

% This function runs the Kalman Filter given parameters. Following Hamilton's (1994) notation
% 
%   State Equation
%   xi_{t+1} = F xi_{t} + v_{t+1}
%
%   Measurement Eq
%   y_t = A' x_t + H' xi_{t} + w_t
%
%   Error covariance
%   Q = E[v_t v_t']
%   R = E[w_t w_t']
%
% INPUTS
% y = T x n_y matrix of n_y measurment time series observations
% x = T x n_x matrix of n_x exogenous variables. This includes constants and lags.
% F = n_s x n_s matrix transition for states 
% A = n_x x n_y matrix for exogenous variables
% H = n_s x n_y matrix for states variables
% Q = n_s x n_s, var-cov matrix for innovations in the state eq.
% R = n_y x n_y, var-cov matrix for innovations in the measurement eq.
% settings: object that contains several options



[T_y,n_y] = size(y);
[T_x,n_x] = size(x);
n_s = size(F,1);

if (T_y ~= T_x); error('y and x have different number of observations'); end;

T = T_y;

% create placeholders for the inferred variables

% xi_{t|t}, contemporanous prediction
xi_mat_c = zeros(T,n_s);
% xi_{t+1|t}, one step ahead prediction
xi_mat = zeros(T,n_s);

% xi_{t|T}, smoothed
xi_smoothed_mat = zeros(T,n_s);

% P_{t|t}
P_mat_c = zeros(T,n_s,n_s);
% P_{t+1|t}
P_mat = zeros(T,n_s,n_s);
% P_{t|T}, smoothed
P_smoothed_mat = zeros(T,n_s,n_s);

% J_{t}
J_mat = zeros(T,n_s,n_s);


loglik_vec = zeros(T,1);


% initial value of the state
try
    xi_0 = settings.xi_0';
catch
    xi_0 = zeros(n_s,1)';
end

% initial value of the precision
try
   P_0   = settings.P_0;
catch
    P_0  = solve_ricatti(F,Q);
end

for t=1:T
    if t == 1
        P_tt = P_0 ;
        xi_tt = xi_0';
    else
        P_tt  = squeeze(P_mat_c(t-1,:,:));
        xi_tt = xi_mat_c(t-1,:)';
    end
    R = R_base*(kappa_seq(t)^2); %this accounts for the time-varying matrix

    %update state one step ahead forecast
    xi_mat(t,:) = (F*xi_tt)';

    %update one step ahead forecast prediction error
    P_ttm1 =  F*P_tt*F'+Q;
    P_mat(t,:,:) = P_ttm1;
    
    % Prediction error
    pred_error = y(t,:)'-(A')*(x(t,:)')-(H')*(xi_mat(t,:)');
    HPHR = (H')*P_ttm1 *H+R;

    % Given this compute the likelihood of the obs using eq 13.4.1 in Hamilton
    loglik_vec(t) = -(n_y/2)*log(2*pi)-0.5*log(det(HPHR))...
        -0.5*pred_error'*(HPHR\pred_error);

    % Contemporaneous prediction
    xi_mat_c(t,:) = xi_mat(t,:)'+...
        (P_ttm1*H)*(HPHR\pred_error);

     % Contemporaneous prediction error:
    P_mat_c(t,:,:) = P_ttm1-P_ttm1*H*(HPHR\(H'*P_ttm1));
end

% given this, now run the Kalman smoother to obtain smoothed state
% estimates.

% CHECK THIS AGAINST HLW

if settings.smoother == 1

states_smoother = settings.states_smoother;
n_smoother = length(states_smoother);
% start the smoother
xi_smoothed_mat(T,states_smoother) =  xi_mat_c(T,states_smoother);
P_smoothed_mat(T,states_smoother,states_smoother) = squeeze(P_mat_c(T,states_smoother,states_smoother));
for t = T-1:-1:1
    P_tt = squeeze(P_mat_c(t,states_smoother,states_smoother));
    P_tp1t = squeeze(P_mat(t+1,states_smoother,states_smoother)); 
    %J_t =  P_tt*F(states_smoother,states_smoother)'*(P_tp1t\eye(n_smoother));
    J_t =  P_tt*F(states_smoother,states_smoother)'*pinv(P_tp1t);
    J_mat(t,states_smoother,states_smoother) = J_t;
    xi_tt = (xi_mat_c(t,states_smoother)');
    xt_tp1t = xi_mat(t+1,states_smoother)';
    xi_tT = xi_tt + J_t*((xi_smoothed_mat(t+1,states_smoother)')-xt_tp1t);
    P_tT = P_tt + J_t*(squeeze(P_smoothed_mat(t+1,states_smoother,states_smoother))-P_tp1t)*J_t';
    xi_smoothed_mat(t,states_smoother) = xi_tT;
    P_smoothed_mat(t,states_smoother,states_smoother) = P_tT;
end

end


output.state_contemporaneous = xi_mat_c;
output.state_contemporaneous_pred = P_mat_c;
output.state_one_ahead = xi_mat;
output.state_one_ahead_pred = P_mat;
output.state_smoothed = xi_smoothed_mat;
output.state_smoothed_pred = P_smoothed_mat;
output.loglik_vec = loglik_vec;
