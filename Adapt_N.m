function [N_next,e_prev] = Adapt_N(feedback, prediction, e_curr, N_curr, err_tol,err_der_tol,N_min,N_max,Ts)

%%%%%%%%%%%%%%%%%%%%%%%% Prediction Horizon Estimator%%%%%%%%%%%%%%%%%%%

% Compute the error and its derivative
e            = prediction - feedback;
de           = (e - e_curr)/Ts;
e_prev       = e;

% Adapt prediction horizon only based on the "Position"

if norm(e(1)) > err_tol || norm(de(1)) > err_der_tol % Reduce prediction horizon if error is large
    N_next        = max(N_min, N_curr - 1);
else % Increase prediction horizon if error is small
    N_next         = min(N_max, N_curr + 1);


end