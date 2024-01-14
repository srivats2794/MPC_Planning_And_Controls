function [pl_solver,pl_args,f] = pl_prob_setup_kin(pl,sim,sys)
% Function that sets up the planning problem

import casadi.*;

x = SX.sym('x'); y = SX.sym('y'); psi = SX.sym('psi');

states = [x;y;psi]; n_states = length(states);

v_l = SX.sym('v_l'); v_r = SX.sym('v_r');
controls = [v_l;v_r]; n_controls = length(controls);

rhs = [ ((v_l+v_r)/2)*cos(psi); ...
    ((v_l+v_r)/2)*sin(psi); ...
    ((v_l-v_r)/sys.w);];

f = Function('f',{states,controls},{rhs});
U = SX.sym('U',n_controls,pl.N);
P = SX.sym('P',n_states + pl.nxo +pl.nu);
X = SX.sym('X',n_states,(pl.N+1));

obj = 0; % Objective function
g = [];  % constraints vector
g_in=[];
st  = X(:,1); % initial state
g = [g;st-P(1:n_states)]; % initial condition constraints
con_ref= P(7:8);
for k = 1:pl.N
    st = X(:,k);  con = U(:,k);

    obj_term1= (st(1:pl.nxo)-P(n_states+1:n_states+pl.nxo));
    obj_term3= con;


    % Calculate objective

    obj = obj+obj_term1'*pl.Q*obj_term1+ ...
             +obj_term3'*pl.R*obj_term3;

    % RK4 Integration Scheme
    st_next = X(:,k+1);
    k1 = f(st, con);
    k2 = f(st + (pl.dt/2)*k1, con);
    k3 = f(st + (pl.dt/2)*k2, con);
    k4 = f(st + pl.dt*k3, con);
    st_Dot = (1/6)*(k1 +2*k2 +2*k3 +k4);
    st_next_RK4=st +pl.dt*st_Dot;
    %

    g = [g;st_next-st_next_RK4]; % Compute equality constraints

    g_in= [g_in;(con-con_ref)/pl.dt;st_Dot(3)];
end
g=[g;g_in];
obj_term_terminal= (st(1:2)-P(n_states+1:n_states+2)); % Terminal Constraint
obj= obj+obj_term_terminal'*pl.QE*obj_term_terminal;

% Add constraints for collision avoidance
for i=1:sim.obs_num
    for k = 1:pl.N+1
        g = [g ; -sqrt((X(1,k)-sim.obs_x(i))^2+(X(2,k)-sim.obs_y(i))^2) + pl.ego_safety_radius/2 + sim.obs_diam/2];
    end
end

% make the decision variable one column  vector
OPT_variables = [reshape(X,n_states*(pl.N+1),1);reshape(U,n_controls*pl.N,1)];

nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);

opts = struct;
opts.ipopt.max_iter = 100;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-2;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

pl_solver = nlpsol('solver', 'ipopt', nlp_prob,opts);

pl_args = struct;
pl_args.lbg(1:n_states*(pl.N+1)) = 0; % equality constraints
pl_args.ubg(1:n_states*(pl.N+1)) = 0; % equality constraints
pl_args.lbg(n_states*(pl.N+1)+1 : n_states: n_states*(pl.N+1) + n_states*pl.N) = pl.a_min; % inequality constraints
pl_args.lbg(n_states*(pl.N+1)+2 : n_states: n_states*(pl.N+1) + n_states*pl.N) = pl.a_min;   % inequality constraints
pl_args.ubg(n_states*(pl.N+1)+1 : n_states: n_states*(pl.N+1) + n_states*pl.N) = pl.a_max; % inequality constraints
pl_args.ubg(n_states*(pl.N+1)+2 : n_states: n_states*(pl.N+1) + n_states*pl.N) = pl.a_max;   % inequality constraints
pl_args.lbg(n_states*(pl.N+1)+3 : n_states: n_states*(pl.N+1) + n_states*pl.N) = pl.psi_dot_min; % inequality constraints
pl_args.ubg(n_states*(pl.N+1)+3 : n_states: n_states*(pl.N+1) + n_states*pl.N) = pl.psi_dot_max; % inequality constraints

pl_args.lbg(n_states*(pl.N+1) + n_states*pl.N+1 : n_states*(pl.N+1) + n_states*pl.N+ sim.obs_num*(pl.N+1)) = -inf; % inequality constraints
pl_args.ubg(n_states*(pl.N+1) + n_states*pl.N+1 : n_states*(pl.N+1) + n_states*pl.N+ sim.obs_num*(pl.N+1)) = 0; % inequality constraints

pl_args.lbx(1:n_states:n_states*(pl.N+1),1) = 0.25; %state x lower bound
pl_args.ubx(1:n_states:n_states*(pl.N+1),1) = 2.75; %state x upper bound
pl_args.lbx(2:n_states:n_states*(pl.N+1),1) = 0; %state y lower bound
pl_args.ubx(2:n_states:n_states*(pl.N+1),1) = 4.25; %state y upper bound
pl_args.lbx(3:n_states:n_states*(pl.N+1),1) = -inf; %state theta lower bound
pl_args.ubx(3:n_states:n_states*(pl.N+1),1) = inf; %state theta upper bound

pl_args.lbx(n_states*(pl.N+1)+1:n_controls:n_states*(pl.N+1)+n_controls*pl.N,1) = pl.v_min; %V_L lower bound
pl_args.ubx(n_states*(pl.N+1)+1:n_controls:n_states*(pl.N+1)+n_controls*pl.N,1) = pl.v_max; %V_L upper bound
pl_args.lbx(n_states*(pl.N+1)+2:n_controls:n_states*(pl.N+1)+n_controls*pl.N,1) = pl.v_min; %V_R lower bound
pl_args.ubx(n_states*(pl.N+1)+2:n_controls:n_states*(pl.N+1)+n_controls*pl.N,1) = pl.v_max; %V_R upper bound
end

