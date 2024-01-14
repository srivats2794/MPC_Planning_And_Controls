clc;clear all;

%% System Parameters
sys.w = 0.44;                               % Width
sys.r_w = 0.1;                              % Radius of Wheel
sys.j_m_th=0.153;                           % Base pitch moment of inertia
sys.j_p_th= 0.125;                          % Pendulum pitch moment of interia    
sys.j_w= 0.013;                             % Rotational inertia of wheel
sys.m_b= 15.747;                            % Base mass
sys.m_w= 5.44;                              % Wheel mass
sys.m_p=4;                                  % Pendulum mass
sys.l= 0.53;                                % Pendulum length
sys.j_psi= 0.576;                           % Base yaw moment of inertia
sys.g= 9.8056;

%% Simulation Parameters 
sim.x0= [0.75 ; 0.65 ;0; 0; 0;0;0];
sim.xf= [1.8; 3.25 ; pi/2; 0; 0;0;0];
sim.obs_num=4; sim.obs_diam = 0.6; 
sim.obs_x= [0.9;2;1.4;2.5;];
sim.obs_y= [1.5;1.75;2.8;3];
% sim.obs_x= [0.3;0.9;1.5;2.1;];
% sim.obs_y= [1.5;1.5;1.5;1.5];
sim.tsim=  10;
sim.x_min= 0;
sim.x_max= 3;
sim.y_min= 0;
sim.y_max= 5;

%% Planner Initialization - NMPC Planner
pl.choice = 2; % 1 Kinematic Diff Drive, 2 Dynamic Diff Drive, 3 Dynamic Self Balancing

pl.Ts = 0.1; % Sampling Rate 
pl.dt = 0.1; % Prediction Interval

if pl.choice == 1
    pl.Q  = diag([0.01;0.01;0.002]); % Penalty Position Regular
    pl.N  = 20; % Prediction Horizon
    pl.R  = diag([0.1;0.1]); % Penalty Input
    pl.QE= diag([1;1]); % Penalty Position Terminal
    pl.v_max = 4.5; % Max Forward Vel
    pl.v_min = -4.5; % Max Backward Vel
    pl.psi_dot_max = 0.785398; % Max AntiClockwise Yawrate
    pl.psi_dot_min =-0.785398; % Max Clockwise Yawrate
    pl.a_max = 2.9; % Max Accel
    pl.a_min = -2.9; % Max Decel
    pl.ego_safety_radius= 0.5;
    
    pl.nx= 3; % Total number of states
    pl.nu= 2; % Total number of inputs
    pl.nxo=3; % Number of states in objective function
    tic
    [pl_solver,pl_args,f_temp]= pl_prob_setup_kin(pl,sim,sys); % CasADi solver setup
    toc
elseif pl.choice==2
    pl.Q  = diag([0.01;0.01;0.002]); % Penalty Position Regular
    pl.N  = 150; % Prediction Horizon
    pl.R  = diag([0.1;0.1]); % Penalty Input
    pl.QE= diag([1;1]); % Penalty Position Terminal
    pl.pl_sys= ctrl_sys_setup_mpc(sys);
    pl.v_max = 4.5; % Max Forward Vel
    pl.v_min = -4.5; % Max Backward Vel
    pl.psi_dot_max = 0.785398; % Max AntiClockwise Yawrate
    pl.psi_dot_min =-0.785398; % Max Clockwise Yawrate
    pl.a_max = 2.9; % Max Accel
    pl.a_min = -2.9; % Max Decel
    pl.ego_safety_radius= 0.6;
    pl.tau_min= -11.5; % Max Reverse Torque
    pl.tau_max= -pl.tau_min; % Max Forward Torque
    pl.nx= 5; % Total number of states
    pl.nu= 2; % Total number of inputs
    pl.nxo=3; % Number of states in objective function
    [pl_solver,pl_args,f_temp]= pl_prob_setup_kinodyn(pl,sim,sys); % CasADi solver setup
elseif pl.choice == 3
    pl.Q  = diag([0.1;0.1]); % Penalty Position Regular
    pl.N  = 200; % Prediction Horizon
    pl.R  = diag([0.1;0.1]); % Penalty Input
    pl.Q2 = diag([0.1;0.05]); % Penalty Pitch
    pl.QE= diag([0;0]); % Penalty Position Terminal
    pl.v_max = 4.5; % Max Forward Vel
    pl.v_min = -4.5; % Max Backward Vel
    pl.psi_dot_max = 0.785398; % Max AntiClockwise Yawrate
    pl.psi_dot_min =-0.785398; % Max Clockwise Yawrate
    pl.a_max = 2.9; % Max Accel
    pl.a_min = -2.9; % Max Decel
    pl.ego_safety_radius= 0.5;
    pl.tau_min= -11.5; % Max Reverse Torque
    pl.tau_max= -pl.tau_min; % Max Forward Torque

    pl.nx= 7; % Total number of states
    pl.nu= 2; % Total number of inputs
    pl.nxo=2; % Number of states in objective function
    tic
    [pl_solver,pl_args,f_temp]= pl_prob_setup(pl,sim,sys); % CasADi solver setup
    toc
end
%% Controller Initialization
if pl.choice==1
    ctrl.Ts= 0.01;  % Sampling Rate
    ctrl.tau_min= -11.5; % Max Reverse Torque
    ctrl.tau_max= -ctrl.tau_min; % Max Forward Torque
    ctrl.ctrl_sys= ctrl_sys_setup_mpc(sys);

    ctrl.nx= length(ctrl.ctrl_sys.A);
    [~,ctrl.nu]= size(ctrl.ctrl_sys.B);
    sys_c= ss(ctrl.ctrl_sys.A,ctrl.ctrl_sys.B,ctrl.ctrl_sys.C,ctrl.ctrl_sys.D);
    [ctrl.ctrl_sys.Ad,ctrl.ctrl_sys.Bd,~,~]=ssdata(c2d(sys_c,ctrl.Ts));
    clearvars sys_c
    ctrl.Q= diag([1;1;2;2]);
    ctrl.R= diag([5;5]);
    ctrl.N=100;
    ctrl.x_max= [pl.v_max;pl.v_max;0.1;0.1];
    ctrl.x_min= [pl.v_min;pl.v_min;-0.1;-0.1];
    ctrl.solver= ctrl_prob_setup_mpc(ctrl.ctrl_sys,ctrl);
elseif pl.choice == 2
    ctrl.Ts= 0.005;
    ctrl.ctrl_sys= ctrl_sys_setup_mpc(sys);
    sys_c= ss(ctrl.ctrl_sys.A_th,ctrl.ctrl_sys.B_th,ctrl.ctrl_sys.C_th,ctrl.ctrl_sys.D_th);
    [ctrl.ctrl_sys.Ad_th,ctrl.ctrl_sys.Bd_th,~,~]=ssdata(c2d(sys_c,ctrl.Ts));
    clearvars sys_c
    ctrl.Q= diag([0.1;1]);
    ctrl.R= diag([1;1]);
    [~,ctrl.K,~] = idare(ctrl.ctrl_sys.Ad_th,ctrl.ctrl_sys.Bd_th,ctrl.Q,ctrl.R,[],[]);
    ctrl.K= -ctrl.K;
elseif pl.choice == 3
end
%% Viz Setup - For simulation graphics visualization
viz.w= sys.w;
viz.l= sys.r_w*2+0.02;
viz.marker_h= viz.l;
viz.marker_w= viz.w/4;
viz.map_x= [0, 3, 3, 0, 0];
viz.map_y= [0, 0, 4.5, 4.5, 0];
viz.ang= 0:0.005:2*pi;
viz.obs_circ_x= sim.obs_diam/2*cos(viz.ang);
viz.obs_circ_y= sim.obs_diam/2*sin(viz.ang);
viz.rob_circ_x= (sys.w/2+0.01)*cos(viz.ang);
viz.rob_circ_y= (sys.w/2+0.01)*sin(viz.ang);

%% Run Simulation
if pl.choice==1
    main_nlmpc;
elseif pl.choice == 2
    main_nmpc_lqr;
elseif pl.choice == 3
    main;
end
%% Run Visualization
sim_viz(pl,sim,viz,pos_fbk_vec,pl_rec);