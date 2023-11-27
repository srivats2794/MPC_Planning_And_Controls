t0 = 0;

pos_fbk_vec(:,1) = sim.x0; % xx contains the history of states
fbk=sim.x0;
t(1) = t0;

u0_pl = zeros(pl.N,2);        % two control inputs for each robot
X0_pl = repmat(sim.x0,1,pl.N+1)'; % initialization of the states decision variables

u_ref=[0;0];
% Start MPC
pl_rec = [];
ctrl_ref=[];


main_loop = tic;
dt_curr= pl.dt;
K= floor(sim.tsim/ctrl.Ts);
pl_exec_freq= round(pl.Ts/ctrl.Ts);
for i=1:K
    if(norm((fbk(1:3)-sim.xf),2) < 5e-3)
        break;
    end

    if(rem(i,pl_exec_freq)==0) %% Planning Cycle
       
        pl_args.p   = [fbk(1:3);sim.xf;u_ref];
        pl_args.x0  = [reshape(X0_pl',3*(pl.N+1),1);reshape(u0_pl',2*pl.N,1)];
        pl_sol = pl_solver('x0', pl_args.x0, 'lbx', pl_args.lbx, 'ubx', pl_args.ubx,...
            'lbg', pl_args.lbg, 'ubg', pl_args.ubg,'p',pl_args.p);
        pl_u = reshape(full(pl_sol.x(3*(pl.N+1)+1:end))',2,pl.N)'; % get controls only from the solution
        pl_rec(:,1:3,i)= reshape(full(pl_sol.x(1:3*(pl.N+1)))',3,pl.N+1)'; % get solution TRAJECTORY
        ctrl_ref= [ctrl_ref ; pl_u(1,:)];
        t(i) = t0;
        % Apply the control and shift the solution
        [t0, fbk, u0_pl] = shift(pl.Ts, t0, fbk, pl_u,f_temp);
        pos_fbk_vec(:,i+1) = fbk(1:3);
        X0_pl = reshape(full(pl_sol.x(1:3*(pl.N+1)))',3,pl.N+1)'; % get solution TRAJECTORY
        % Shift trajectory to initialize the next step
        X0_pl = [X0_pl(2:end,:);X0_pl(end,:)];
        u_ref=pl_u(1,:)';
    end

    
end
main_loop_time = toc(main_loop);


ss_error = norm((fbk-sim.xf),2)

average_mpc_time = main_loop_time/(i)

% t_ref= (0:0.1:24.9)';
% t_vec= (0:0.01:24.99)';
% 
% planner_out_fit(:,1)= makima(t_ref,planner_out(:,1),t_vec);
% planner_out_fit(:,2)= makima(t_ref,planner_out(:,2),t_vec);
% 
% nexttile
% plot(t_ref,planner_out(:,1));
% hold on
% plot(t_vec,planner_out_fit(:,1));
% hold off
% nexttile
% plot(t_ref,planner_out(:,2));
% hold on
% plot(t_vec,planner_out_fit(:,2));
% hold off
% nexttile
% plot(t_ref(2:end),diff(planner_out(:,1))./diff(t_ref));
% % [n,~]=size(planner_out);
% % Tvec= 0:pl.Ts:((n-1)*pl.Ts);
% % figure(2)
% % plot(Tvec,planner_out(:,1),'-k',Tvec,planner_out(:,2),'--r','LineWidth',2);
% % legend('Left Wheel Velocity (m/s)','Right Wheel Velocity (m/s)','FontSize',10);
