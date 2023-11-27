function [pl_solver,pl_args,f] = pl_prob_setup(pl,sim,sys)
import casadi.*;
x = SX.sym('x'); y = SX.sym('y'); psi = SX.sym('psi');
states = [x;y;psi]; n_states = length(states);

V_L = SX.sym('V_L'); V_R = SX.sym('V_R');
controls = [V_L;V_R]; n_controls = length(controls);
rhs = [((V_L+V_R)/2)*cos(psi);((V_L+V_R)/2)*sin(psi);((V_L-V_R)/sys.w)]; 

f = Function('f',{states,controls},{rhs}); 
U = SX.sym('U',n_controls,pl.N); 
P = SX.sym('P',n_states + n_states + n_controls);
X = SX.sym('X',n_states,(pl.N+1));

obj = 0; % Objective function
g = [];  % constraints vector

st  = X(:,1); % initial state
g = [g;st-P(1:3)]; % initial condition constraints
g_in= [];
con_ref= P(7:8);
for k = 1:pl.N
    st = X(:,k);  con = U(:,k);
    obj = obj+(st-P(4:6))'*pl.Q*(st-P(4:6))+con'*pl.R*con; % calculate obj
    st_next = X(:,k+1);
    k1 = f(st, con);   % new 
    k2 = f(st + pl.dt/2*k1, con); % new
    k3 = f(st + pl.dt/2*k2, con); % new
    k4 = f(st + pl.dt*k3, con); % new
    st_next_RK4=st +pl.dt/6*(k1 +2*k2 +2*k3 +k4); % new    
    
    g = [g;st_next-st_next_RK4]; % compute constraints % new
    g_in= [g_in;(con-con_ref)/pl.dt;k1(3)];
    U_ref= con;
end

g=[g;g_in];

% Add constraints for collision avoidance
for i=1:sim.obs_num
    for k = 1:pl.N+1   % box constraints due to the map margins
        g = [g ; -sqrt((X(1,k)-sim.obs_x(i))^2+(X(2,k)-sim.obs_y(i))^2) + pl.ego_safety_radius/2 + sim.obs_diam/2];
    end
end

% make the decision variable one column  vector
OPT_variables = [reshape(X,3*(pl.N+1),1);reshape(U,2*pl.N,1)];

nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);

opts = struct;
opts.ipopt.max_iter = 100;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

pl_solver = nlpsol('solver', 'ipopt', nlp_prob,opts);

pl_args = struct;
pl_args.lbg(1:3*(pl.N+1)) = 0; % equality constraints
pl_args.ubg(1:3*(pl.N+1)) = 0; % equality constraints

pl_args.lbg(3*(pl.N+1)+1 : 3: 3*(pl.N+1) + 3*pl.N) = pl.a_min; % inequality constraints
pl_args.lbg(3*(pl.N+1)+2 : 3: 3*(pl.N+1) + 3*pl.N) = pl.a_min;   % inequality constraints
pl_args.ubg(3*(pl.N+1)+1 : 3: 3*(pl.N+1) + 3*pl.N) = pl.a_max; % inequality constraints
pl_args.ubg(3*(pl.N+1)+2 : 3: 3*(pl.N+1) + 3*pl.N) = pl.a_max;   % inequality constraints
pl_args.lbg(3*(pl.N+1)+3 : 3: 3*(pl.N+1) + 3*pl.N) = pl.psi_dot_min; % inequality constraints
pl_args.ubg(3*(pl.N+1)+3 : 3: 3*(pl.N+1) + 3*pl.N) = pl.psi_dot_max; % inequality constraints

pl_args.lbg(3*(pl.N+1) + 3*pl.N+1 : 3*(pl.N+1) + 3*pl.N+ sim.obs_num*(pl.N+1)) = -inf; % inequality constraints
pl_args.ubg(3*(pl.N+1) + 3*pl.N+1 : 3*(pl.N+1) + 3*pl.N+ sim.obs_num*(pl.N+1)) = 0; % inequality constraints

pl_args.lbx(1:3:3*(pl.N+1),1) = sim.x_min; %state x lower bound
pl_args.ubx(1:3:3*(pl.N+1),1) = sim.x_max; %state x upper bound
pl_args.lbx(2:3:3*(pl.N+1),1) = sim.y_min; %state y lower bound
pl_args.ubx(2:3:3*(pl.N+1),1) = sim.y_max; %state y upper bound
pl_args.lbx(3:3:3*(pl.N+1),1) = -inf; %state theta lower bound
pl_args.ubx(3:3:3*(pl.N+1),1) = inf; %state theta upper bound

pl_args.lbx(3*(pl.N+1)+1:2:3*(pl.N+1)+2*pl.N,1) = pl.v_min; %V_L lower bound
pl_args.ubx(3*(pl.N+1)+1:2:3*(pl.N+1)+2*pl.N,1) = pl.v_max; %V_L upper bound
pl_args.lbx(3*(pl.N+1)+2:2:3*(pl.N+1)+2*pl.N,1) = pl.v_min; %V_R lower bound
pl_args.ubx(3*(pl.N+1)+2:2:3*(pl.N+1)+2*pl.N,1) = pl.v_max; %V_R upper bound
end

