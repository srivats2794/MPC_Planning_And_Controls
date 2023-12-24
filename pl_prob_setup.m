function [pl_solver,pl_args,f] = pl_prob_setup(pl,sim,sys)
import casadi.*;

%% Predefining some repeated cluster terms
    Me= sys.m_b+2*sys.m_w+sys.m_p+((2*sys.j_w)/(sys.r_w^2));
    Je= sys.j_m_th+sys.j_p_th;
    
    %% Predifining some big constants to avoid typos
    p01= sys.m_p*sys.l*sys.l+Je;
    p04= sys.m_p*sys.m_p*sys.l*sys.l;

    p4= p01/(Me*sys.r_w*p01-p04*sys.r_w);
    p5= (sys.w*sys.w*sys.r_w)/ ...
        ((2*sys.j_psi*sys.r_w*sys.r_w)+ ...
        ((sys.j_w+sys.r_w*sys.r_w*sys.m_w)*sys.w*sys.w));

x = SX.sym('x'); y = SX.sym('y'); psi = SX.sym('psi'); v= SX.sym('v');
psiDot= SX.sym('psiDot');
states = [x;y;psi;v;psiDot]; n_states = length(states);

tau_L = SX.sym('tau_L'); tau_R = SX.sym('tau_R');
controls = [tau_L;tau_R]; n_controls = length(controls);
rhs = [v*cos(psi);v*sin(psi);psiDot;p4*(tau_L+tau_R);(2*p5/sys.w)*(tau_L-tau_R)]; 

f = Function('f',{states,controls},{rhs}); 
U = SX.sym('U',n_controls,pl.N); 
P = SX.sym('P',n_states + n_states-2);
X = SX.sym('X',n_states,(pl.N+1));

obj = 0; % Objective function
g = [];  % constraints vector

st  = X(:,1); % initial state
g = [g;st-P(1:5)]; % initial condition constraints

for k = 1:pl.N
    st = X(:,k);  con = U(:,k);
    obj = obj+(st(1:3)-P(6:8))'*pl.Q*(st(1:3)-P(6:8))+con'*pl.R*con; % calculate obj
    st_next = X(:,k+1);
    k1 = f(st, con);   % new 
    k2 = f(st + pl.dt/2*k1, con); % new
    k3 = f(st + pl.dt/2*k2, con); % new
    k4 = f(st + pl.dt*k3, con); % new
    st_next_RK4=st +pl.dt/6*(k1 +2*k2 +2*k3 +k4); % new    
    
    g = [g;st_next-st_next_RK4]; % compute constraints % new
end

% Add constraints for collision avoidance
for i=1:sim.obs_num
    for k = 1:pl.N+1   % box constraints due to the map margins
        g = [g ; -sqrt((X(1,k)-sim.obs_x(i))^2+(X(2,k)-sim.obs_y(i))^2) + pl.ego_safety_radius/2 + sim.obs_diam/2];
    end
end

% make the decision variable one column  vector
OPT_variables = [reshape(X,5*(pl.N+1),1);reshape(U,2*pl.N,1)];

nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);

opts = struct;
opts.ipopt.max_iter = 100;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

pl_solver = nlpsol('solver', 'ipopt', nlp_prob,opts);

pl_args = struct;
pl_args.lbg(1:5*(pl.N+1)) = 0; % equality constraints
pl_args.ubg(1:5*(pl.N+1)) = 0; % equality constraints

pl_args.lbg(5*(pl.N+1) + 1 : length(g)) = -inf; % inequality constraints
pl_args.ubg(5*(pl.N+1) + 1 : length(g)) = 0; % inequality constraints

pl_args.lbx(1:5:5*(pl.N+1),1) = sim.x_min; %state x lower bound
pl_args.ubx(1:5:5*(pl.N+1),1) = sim.x_max; %state x upper bound
pl_args.lbx(2:5:5*(pl.N+1),1) = sim.y_min; %state y lower bound
pl_args.ubx(2:5:5*(pl.N+1),1) = sim.y_max; %state y upper bound
pl_args.lbx(3:5:5*(pl.N+1),1) = -inf; %state psi lower bound
pl_args.ubx(3:5:5*(pl.N+1),1) = inf; %state psi upper bound
pl_args.lbx(4:5:5*(pl.N+1),1) = pl.v_min; %state v lower bound
pl_args.ubx(4:5:5*(pl.N+1),1) = pl.v_max; %state v upper bound
pl_args.lbx(5:5:5*(pl.N+1),1) = pl.psi_dot_min; %state psiDot lower bound
pl_args.ubx(5:5:5*(pl.N+1),1) = pl.psi_dot_max; %state psiDot upper bound

pl_args.lbx(5*(pl.N+1)+1:2:5*(pl.N+1)+2*pl.N,1) = pl.tau_min; %V_L lower bound
pl_args.ubx(5*(pl.N+1)+1:2:5*(pl.N+1)+2*pl.N,1) = pl.tau_max; %V_L upper bound
pl_args.lbx(5*(pl.N+1)+2:2:5*(pl.N+1)+2*pl.N,1) = pl.tau_min; %V_R lower bound
pl_args.ubx(5*(pl.N+1)+2:2:5*(pl.N+1)+2*pl.N,1) = pl.tau_max; %V_R upper bound
end

