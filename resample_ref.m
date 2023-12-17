function ref_resampled = resample_ref(N_c,N_p,dt_c,dt_p,ref_og,t_beg)
% Resamples the planner generated reference to the discretization of the controller

t_pl= (linspace(0,(N_p*dt_p)-dt_p,N_p))';
t_ctrl= (linspace(t_beg,t_beg+(N_c*dt_c)-dt_c,N_c))';

ref_resampled = makima(t_pl,ref_og,t_ctrl);

end

