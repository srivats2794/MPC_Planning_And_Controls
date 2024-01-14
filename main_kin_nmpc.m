fbk=sim.x0(1:pl.nx);

% two inputs for robot at planner level
u0_pl = zeros(pl.N,2);       
% initialization of the kinematic states decision variables
X0_pl = repmat(sim.x0(1:pl.nx),1,pl.N+1)'; 

u_ref=[0;0];

pl_rec = [];
ctrl_ref=[];

K= floor(sim.tsim/pl.Ts);
main_loop = tic;
for i=1:K
    if(norm((fbk(1:pl.nxo)-sim.xf(1:pl.nxo)),2) < 0.1)
        break;
    end
    norm((fbk(1:pl.nxo)-sim.xf(1:pl.nxo)),2)
    t(i) = (i-1)*pl.Ts;

    
    pl_args.p   = [fbk;sim.xf(1:pl.nxo);u_ref];
    pl_args.x0  = [reshape(X0_pl',pl.nx*(pl.N+1),1);reshape(u0_pl',pl.nu*pl.N,1)];
    pl_sol = pl_solver('x0', pl_args.x0, 'lbx', pl_args.lbx, 'ubx', pl_args.ubx,...
        'lbg', pl_args.lbg, 'ubg', pl_args.ubg,'p',pl_args.p);

    % get inputs only from the solution
    pl_u = reshape(full(pl_sol.x(pl.nx*(pl.N+1)+1:end))',pl.nu,pl.N)';
    % get solution TRAJECTORY for plotting purposes
    pl_st= reshape(full(pl_sol.x(1:pl.nx*(pl.N+1)))',pl.nx,pl.N+1)';
    pl_rec(:,1:pl.nx,i)= pl_st(:,1:3);

    X0_pl = reshape(full(pl_sol.x(1:pl.nx*(pl.N+1)))',pl.nx,pl.N+1)'; % get solution TRAJECTORY

    % Shift trajectory to initialize the next step
    X0_pl = [X0_pl(2:end,:);X0_pl(end,:)];
    rel_t = pl.dt:pl.dt:pl.dt*(pl.N+1);

    act_t= pl.Ts:pl.dt:pl.dt*(pl.N+1);

    X0_pl= makima(rel_t',X0_pl',act_t');
    X0_pl= X0_pl';
    % record fbk for plotting purposes
    pos_fbk_vec(:,i) = fbk(1:3);

    % dyn_fbk_vec(:,i) = [fbk(4); ...
    %     fbk(5);
    %     fbk(6);
    %     fbk(7);];

    % tau_l(i)= pl_u(1,1);
    % tau_r(i)= pl_u(1,2);
    v_l(i)= pl_u(1,1);
    v_r(i)= pl_u(1,2);

    %fbk= propagate_plant(sys,fbk,[tau_l(i);tau_r(i)],pl.Ts,0,ctrl_sys_setup_mpc(sys));

    k1 = f_temp(fbk, [v_l(i);v_r(i)]);
    k2 = f_temp(fbk + (pl.Ts/2)*k1, [v_l(i);v_r(i)]);
    k3 = f_temp(fbk + (pl.Ts/2)*k2, [v_l(i);v_r(i)]);
    k4 = f_temp(fbk + pl.Ts*k3, [v_l(i);v_r(i)]);
    fbk=full(fbk +(pl.Ts/6)*(k1 +2*k2 +2*k3 +k4));
    
    u_ref= [v_l(i);v_r(i)];
end
main_loop_time = toc(main_loop);

ss_error = norm((fbk(1:2)-sim.xf(1:2)),2)

average_mpc_time = main_loop_time/(i)