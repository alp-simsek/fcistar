function  se_output = compute_se_kalman(theta_opt,constraints,evaluate_kalman,N_draws)

% This script computes the standard error of parameters and states
% (accounting for parameter estimation error) using the OPG to get the
% variance-covariance matrix of the MLE estimator, and then using the
% monte carlo procedure in Hamilton (1986) and also used by HLW to compute
% the SEs for the states.

% INPUTS

% theta_opt N x 1 vector of the MLE estimate
% constraints N x 1 vector of logical values. If an element is = 1 it means
% the parameter hit the constraint and thus we don't report SEs for it.
% evaluate_kalman: function that computes the Kalman Filter. We are
% interested in the log-likelihood vector. This should be only a function
% of theta.

% First, compute the variance covariance matrix using OPG.
UB = constraints.UB;
LB = constraints.LB;
try
    const_f = constraints.const_f;
catch
    const_f = @(x) -10;
end

cons_vec = (theta_opt==LB)+(theta_opt==UB);

unconst_cols = find(cons_vec==0);
const_cols = find(cons_vec==1);
N_param = length(theta_opt);
out_base = evaluate_kalman(theta_opt);
[T,N_states] = size(out_base.state_contemporaneous); 

%variance covariance matrix ONLY for the unconstrained parameters.
V_mat = compute_opg(theta_opt,constraints,evaluate_kalman);

C_mat_aux = chol(V_mat,'lower'); %this satisfies C_mat*C_mat*; 
C_mat = zeros(N_param,N_param);
C_mat(unconst_cols,unconst_cols) = C_mat_aux;
se_vec = NaN(N_param,1);
se_vec(unconst_cols) = sqrt(diag(V_mat));
t_stat_vec =  NaN(N_param,1);
t_stat_vec(unconst_cols) = theta_opt(unconst_cols)./se_vec(unconst_cols);

% compute SEs for the states using the estimator in Hamilton 1986.

% generate parameter draws
theta_draws = theta_opt + C_mat*randn(N_param,N_draws);

% if sum(cons_vec)>0
%     theta_draws = NaN(N_param,N_draws);
%     theta_draws(unconst_cols,:) = theta_draws_aux;
%     theta_draws(const_cols,:) = theta_opt(const_cols); %fix this values at the bound
% else
%     theta_draws = theta_draws_aux;
% end

filter_unc_mat = zeros(T,N_states,N_states);
param_unc_mat = zeros(T,N_states,N_states);

ii=1;
used_draws = 0;
while ii<=N_draws
    theta_i = theta_draws(:,ii);
    % first, check if the draw satisfies the constraints
    if  sum([theta_i>UB; theta_i<LB; (const_f(theta_i)>0)])>0 %if not just go to the next draw
        ii = ii+1;
    else % if it satisfies the constraint, solve the optimizer.
        out_i = evaluate_kalman(theta_i);
        % collect the results. 
        % uncertainty coming from the filter. Smoothed is the only one that
        % makes sense since the MLE uses full sample
        filter_unc_mat = filter_unc_mat + out_i.state_smoothed_pred; 
        % squeeze(out_i.state_smoothed_pred(end,:,:))
        
        % difference of the state from baseline
        state_dif = (out_i.state_smoothed - out_base.state_smoothed);
        % squeeze(out_i.state_smoothed_pred(end,:,:))
        for t=1:T
            param_unc_mat(t,:,:) = squeeze(param_unc_mat(t,:,:))+state_dif(t,:)'*state_dif(t,:); 
        end
         % squeeze(param_unc_mat(end,:,:))

        % once this is done, advance iters.
        used_draws = used_draws+1;
        ii = ii+1;
    end
end

filter_unc = filter_unc_mat/used_draws;
param_unc = param_unc_mat/used_draws;
state_unc_mat = out_base.state_smoothed_pred + param_unc;
%state_unc_mat = filter_unc + param_unc;
%algo = squeeze(filter_unc_mat(end,:,:))
%algo = squeeze(param_unc(end,:,:))
%algo_2 = squeeze(out_base.state_smoothed_pred(end,:,:))
%  squeeze(state_unc_mat(end,:,:))
state_ses = NaN(T,N_states);
state_ses_no_param = NaN(T,N_states);
state_ses_contemp = NaN(T,N_states);
for t=1:T
    state_ses(t,:) = sqrt(diag(squeeze(state_unc_mat(t,:,:))))';
    state_ses_no_param(t,:) = sqrt(diag(squeeze(filter_unc(t,:,:))))';
    state_ses_contemp(t,:) = sqrt(diag(squeeze(out_base.state_contemporaneous_pred(t,:,:))))';
end



se_output.V_mat = V_mat;
se_output.se_vec = se_vec;
se_output.t_stat_vec = t_stat_vec;
se_output.filter_unc = filter_unc; 
se_output.param_unc =  param_unc;
se_output.state_unc_mat =state_unc_mat ;
se_output.state_ses =  state_ses;
se_output.state_ses_no_param =  state_ses_no_param;
se_output.state_ses_contemp =  state_ses_contemp;
end


