fbk=sim.x0;
fbk_vec(:,1) = fbk;

% two inputs for robot at planner level
u0_pl = zeros(pl.N,2);       
% initialization of the kinematic states decision variables
X0_pl = repmat(sim.x0(1:pl.nx),1,pl.N+1)'; 

u_ref=[0;0];

pl_rec = [];
ctrl_ref=[];

K= floor(sim.tsim/ctrl.Ts);

pl_exec_freq= round(pl.Ts/ctrl.Ts);
% Create a relative time vector for planning reference
t_pl= (linspace(0,(pl.N*pl.dt)-pl.dt,pl.N))';
i=1;
for j=1:K
    if(norm((fbk(1:pl.nxo)-sim.xf(1:pl.nxo)),2) < 0.1)
        break;
    end
    norm((fbk(1:pl.nxo)-sim.xf(1:pl.nxo)),2)
    t(j) = (j-1)*ctrl.Ts;

    if(rem(j,pl_exec_freq)==1) %% Planning Cycle
        pl_loop = tic;
        pl_args.p   = [fbk(1:pl.nx);sim.xf(1:pl.nxo);u_ref];
        pl_args.x0  = [reshape(X0_pl',pl.nx*(pl.N+1),1);reshape(u0_pl',pl.nu*pl.N,1)];
        pl_sol = pl_solver('x0', pl_args.x0, 'lbx', pl_args.lbx, 'ubx', pl_args.ubx,...
            'lbg', pl_args.lbg, 'ubg', pl_args.ubg,'p',pl_args.p);

        % get inputs only from the solution
        pl_u = reshape(full(pl_sol.x(pl.nx*(pl.N+1)+1:end))',pl.nu,pl.N)';
        % get solution TRAJECTORY for plotting purposes
        pl_st= reshape(full(pl_sol.x(1:pl.nx*(pl.N+1)))',pl.nx,pl.N+1)';
        pl_rec(:,1:pl.nx,i)= pl_st(:,1:3);

        X0_pl = reshape(full(pl_sol.x(1:pl.nx*(pl.N+1)))',pl.nx,pl.N+1)'; % get solution TRAJECTORY
        
        ctrl_ref= pl_u';

        % Shift trajectory to initialize the next step
        X0_pl = [X0_pl(2:end,:);X0_pl(end,:)];
        % rel_t = pl.dt:pl.dt:pl.dt*(pl.N+1);
        % 
        % act_t= pl.Ts:pl.dt:pl.dt*(pl.N+1);
        % 
        % X0_pl= makima(rel_t',X0_pl',act_t');
        X0_pl= X0_pl';
        
        pl_loop_time(i) = toc(pl_loop);
         % record fbk for plotting purposes
        pos_fbk_vec(:,i) = fbk(1:3);
        v_l(i)= pl_u(1,1);
        v_r(i)= pl_u(1,2);
    
        u_ref= [v_l(i);v_r(i)];
        i=i+1;
        pl_status=1; 
        %pl_ref_curve= makima(t_pl,ctrl_ref); % TODO: Try PCHIP instead
    end
    
    % Whenever planner status goes to 1, our start
    if pl_status==1
        count=0;
    else
        count= count+1;
    end

    % Create a relative time vector with lookahead for the controller
    t_ctrl= (linspace((count+1)*ctrl.Ts,(count+1+ctrl.N)*ctrl.Ts,ctrl.N+1))';

    % Interpolate the planned trajectory to control discretization
    %ctrl_ref_curr= ppval(pl_ref_curve,t_ctrl);
    
    % Reset planner status until next planning execution occurs
    pl_status=0; 
    ctrl_ref_curr= [ctrl_ref(1,2);ctrl_ref(2,2);0;0]*ones(1,ctrl.N+1);
    % theta and thetaDot reference for the Controller MPC is always 0
    
    [tau_l(i),tau_r(i),dyn_fbk_pred(:,i+1)]= controller_mpc(ctrl,ctrl_ref_curr,fbk);

    fbk= propagate_plant(sys,fbk,[tau_l(i);tau_r(i)],pl.Ts,1,ctrl_sys_setup_mpc(sys));
    
end


ss_error = norm((fbk(1:2)-sim.xf(1:2)),2)

average_mpc_time = sum(pl_loop_time)/(i)