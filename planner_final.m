clear all;close all; clc

import casadi.*;

Ts_pl = 0.1; h=0.2;
N = 100; % prediction horizon
rob_diam = 0.44;
x0 = [1.85 ; 0.5 ; 0];    % initial condition.
xs = [1.85; 4.15 ; pi/2]; % Reference posture.
% Add constraints for collision avoidance
obst_num=4; obs_diam = 0.6; 

obs_x= [0.40+0.75;0.30+0.75;0.35+2;0.35+1.5;];
obs_y= [2.2;3.2;2.75;1.25];

Q = zeros(3,3); Q(1,1) =0.1;Q(2,2) = 0.1;Q(3,3) = 0.02; % weighing matrices (states)
R = zeros(2,2); R(1,1) = 2; R(2,2) = 2; % weighing matrices (controls) 

v_max = 4.5; v_min = -v_max; 


x = SX.sym('x'); y = SX.sym('y'); psi = SX.sym('psi');
states = [x;y;psi]; n_states = length(states);

V_L = SX.sym('V_L'); V_R = SX.sym('V_R');
controls = [V_L;V_R]; n_controls = length(controls);
rhs = [((V_L+V_R)/2)*cos(psi);((V_L+V_R)/2)*sin(psi);((V_L-V_R)/rob_diam)]; 

f = Function('f',{states,controls},{rhs}); 
U = SX.sym('U',n_controls,N); 
P = SX.sym('P',n_states + n_states + n_controls);
X = SX.sym('X',n_states,(N+1));

obj = 0; % Objective function
g = [];  % constraints vector

st  = X(:,1); % initial state
g = [g;st-P(1:3)]; % initial condition constraints
g_in= [];
con_ref= P(7:8);
for k = 1:N
    st = X(:,k);  con = U(:,k);
    obj = obj+(st-P(4:6))'*Q*(st-P(4:6))+con'*R*con; % calculate obj
    st_next = X(:,k+1);
    k1 = f(st, con);   % new 
    k2 = f(st + h/2*k1, con); % new
    k3 = f(st + h/2*k2, con); % new
    k4 = f(st + h*k3, con); % new
    st_next_RK4=st +h/6*(k1 +2*k2 +2*k3 +k4); % new    
    
    g = [g;st_next-st_next_RK4]; % compute constraints % new
    g_in= [g_in;(con-con_ref)/h;k1(3)];
    U_ref= con;
end

g=[g;g_in];

% Add constraints for collision avoidance
for i=1:obst_num
    for k = 1:N+1   % box constraints due to the map margins
        g = [g ; -sqrt((X(1,k)-obs_x(i))^2+(X(2,k)-obs_y(i))^2) + (((rob_diam/2)+0.1) + obs_diam/2)];
    end
end

% make the decision variable one column  vector
OPT_variables = [reshape(X,3*(N+1),1);reshape(U,2*N,1)];

nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);

opts = struct;
opts.ipopt.max_iter = 100;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = nlpsol('solver', 'ipopt', nlp_prob,opts);

args = struct;
args.lbg(1:3*(N+1)) = 0; % equality constraints
args.ubg(1:3*(N+1)) = 0; % equality constraints

args.lbg(3*(N+1)+1 : 3: 3*(N+1) + 3*N) = -2.9; % inequality constraints
args.lbg(3*(N+1)+2 : 3: 3*(N+1) + 3*N) = -2.9;   % inequality constraints
args.ubg(3*(N+1)+1 : 3: 3*(N+1) + 3*N) = 2.9; % inequality constraints
args.ubg(3*(N+1)+2 : 3: 3*(N+1) + 3*N) = 2.9;   % inequality constraints
args.lbg(3*(N+1)+3 : 3: 3*(N+1) + 3*N) = -0.5; % inequality constraints
args.ubg(3*(N+1)+3 : 3: 3*(N+1) + 3*N) = 0.5; % inequality constraints

args.lbg(3*(N+1) + 3*N+1 : 3*(N+1) + 3*N+ obst_num*(N+1)) = -inf; % inequality constraints
args.ubg(3*(N+1) + 3*N+1 : 3*(N+1) + 3*N+ obst_num*(N+1)) = 0; % inequality constraints

args.lbx(1:3:3*(N+1),1) = 0.25; %state x lower bound
args.ubx(1:3:3*(N+1),1) = 2.75; %state x upper bound
args.lbx(2:3:3*(N+1),1) = 0; %state y lower bound
args.ubx(2:3:3*(N+1),1) = 4.25; %state y upper bound
args.lbx(3:3:3*(N+1),1) = -inf; %state theta lower bound
args.ubx(3:3:3*(N+1),1) = inf; %state theta upper bound

args.lbx(3*(N+1)+1:2:3*(N+1)+2*N,1) = v_min; %V_L lower bound
args.ubx(3*(N+1)+1:2:3*(N+1)+2*N,1) = v_max; %V_L upper bound
args.lbx(3*(N+1)+2:2:3*(N+1)+2*N,1) = v_min; %V_R lower bound
args.ubx(3*(N+1)+2:2:3*(N+1)+2*N,1) = v_max; %V_R upper bound
%----------------------------------------------
% ALL OF THE ABOVE IS JUST A PROBLEM SET UP


% THE SIMULATION LOOP SHOULD START FROM HERE
%-------------------------------------------
t0 = 0;

planner_hist(:,1) = x0; % xx contains the history of states
t(1) = t0;

u0 = zeros(N,2);        % two control inputs for each robot
X0 = repmat(x0,1,N+1)'; % initialization of the states decision variables

sim_tim = 25; % Maximum simulation time
u_ref=[0;0];
% Start MPC
mpciter = 0;
x_sol_rec = [];
planner_out=[];

% the main simulaton loop... it works as long as the error is greater
% than 10^-6 and the number of mpc steps is less than its maximum
% value.
main_loop = tic;
while(norm((x0-xs),2) > 5e-2 && mpciter < sim_tim / Ts_pl)
    args.p   = [x0;xs;u_ref];
    args.x0  = [reshape(X0',3*(N+1),1);reshape(u0',2*N,1)];
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
        'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
    u = reshape(full(sol.x(3*(N+1)+1:end))',2,N)'; % get controls only from the solution
    x_sol_rec(:,1:3,mpciter+1)= reshape(full(sol.x(1:3*(N+1)))',3,N+1)'; % get solution TRAJECTORY
    planner_out= [planner_out ; u(1,:)];
    t(mpciter+1) = t0;
    % Apply the control and shift the solution
    [t0, x0, u0] = shift(Ts_pl, t0, x0, u,f);
    planner_hist(:,mpciter+2) = x0;
    X0 = reshape(full(sol.x(1:3*(N+1)))',3,N+1)'; % get solution TRAJECTORY
    % Shift trajectory to initialize the next step
    X0 = [X0(2:end,:);X0(end,:)];
    u_ref=u(1,:)';
    mpciter
    mpciter = mpciter + 1;
end;
main_loop_time = toc(main_loop);


ss_error = norm((x0-xs),2)
average_mpc_time = main_loop_time/(mpciter+1)
draw_now;

clearvars -except planner_out Ts_pl

t_ref= (0:0.1:24.9)';
t_vec= (0:0.01:24.99)';

planner_out_fit(:,1)= makima(t_ref,planner_out(:,1),t_vec);
planner_out_fit(:,2)= makima(t_ref,planner_out(:,2),t_vec);

nexttile
plot(t_ref,planner_out(:,1));
hold on
plot(t_vec,planner_out_fit(:,1));
hold off
nexttile
plot(t_ref,planner_out(:,2));
hold on
plot(t_vec,planner_out_fit(:,2));
hold off
nexttile
plot(t_ref(2:end),diff(planner_out(:,1))./diff(t_ref));
% [n,~]=size(planner_out);
% Tvec= 0:Ts_pl:((n-1)*Ts_pl);
% figure(2)
% plot(Tvec,planner_out(:,1),'-k',Tvec,planner_out(:,2),'--r','LineWidth',2);
% legend('Left Wheel Velocity (m/s)','Right Wheel Velocity (m/s)','FontSize',10);
