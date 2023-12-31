fbk=sim.x0;

dyn_fbk_vec(:,1)= [fbk(4)+(sys.w/2)*fbk(5); ...
                   fbk(4)-(sys.w/2)*fbk(5); ...
                   fbk(6); ...
                   fbk(7)];
dyn_pred_vec(:,1)= dyn_fbk_vec(:,1);
% two inputs for robot at planner level
u0_pl = zeros(pl.N,2);       
% initialization of the kinematic states decision variables
X0_pl = repmat(sim.x0(1:3),1,pl.N+1)'; 

u_ref=[0;0];

pl_rec = [];
ctrl_ref=[];

K= floor(sim.tsim/ctrl.Ts);

% Ratio between planning and controls execution cycles
pl_exec_freq= round(pl.Ts/ctrl.Ts);

% Create a relative time vector for planning reference
t_pl= (linspace(0,(pl.N*pl.dt)-pl.dt,pl.N))'; 

main_loop = tic;
j=1;

for i=1:K
    if(norm((fbk(1:3)-sim.xf(1:3)),2) < 5e-3)
        break;
    end
    t(i) = (i-1)*ctrl.Ts;
    if(rem(i,pl_exec_freq)==1) %% Planning Cycle
       
        pl_args.p   = [fbk(1:3);sim.xf(1:3);u_ref];
        pl_args.x0  = [reshape(X0_pl',3*(pl.N+1),1);reshape(u0_pl',2*pl.N,1)];
        pl_sol = pl_solver('x0', pl_args.x0, 'lbx', pl_args.lbx, 'ubx', pl_args.ubx,...
            'lbg', pl_args.lbg, 'ubg', pl_args.ubg,'p',pl_args.p);

        % get inputs only from the solution
        pl_u = reshape(full(pl_sol.x(3*(pl.N+1)+1:end))',2,pl.N)'; 
        % get solution TRAJECTORY for plotting purposes
        pl_rec(:,1:3,i)= reshape(full(pl_sol.x(1:3*(pl.N+1)))',3,pl.N+1)'; 
        % the inputs solved by planner are the references for the
        % controller
        ctrl_ref= (full(pl_u))';
        ctrl_ref= ctrl_ref(:,1);
        
        X0_pl = reshape(full(pl_sol.x(1:3*(pl.N+1)))',3,pl.N+1)'; % get solution TRAJECTORY

        % Shift trajectory to initialize the next step
        X0_pl = [X0_pl(2:end,:);X0_pl(end,:)];

        u_ref=pl_u(1,:)';
        
        % record feedback for plotting purposes
        pos_fbk_vec(:,j) = fbk(1:3);
        j=j+1;
        % A flag that marks that a planner cycle has just gotten done
        pl_status=1; 
        % Store the planner trajectory as a smooth C2 spline to be used by
        % controller
        %pl_ref_curve= makima(t_pl,ctrl_ref); % TODO: Try PCHIP instead
    end
    
    % % Whenever planner status goes to 1, our start
    % if pl_status==1
    %     count=0;
    % else
    %     count= count+1;
    % end
    % 
    % % Create a relative time vector with lookahead for the controller
    % t_ctrl= (linspace((count+ctrl.lookahead)*ctrl.Ts,count*ctrl.Ts+(ctrl.N*ctrl.Ts),ctrl.N+1))';

    % Interpolate the planned trajectory to control discretization
    %ctrl_ref_curr= ppval(pl_ref_curve,t_ctrl);
    
    % Reset planner status until next planning execution occurs;
    %pl_status=0; 
    
    % theta and thetaDot reference for the Controller MPC is always 0
    ctrl_ref(3:4,1)=0;
    
    % Store dynamic states for plotting purposes    
    dyn_fbk_vec(:,i+1) = [fbk(4)+(sys.w/2)*fbk(5); ...
                          fbk(4)-(sys.w/2)*fbk(5); ...
                          fbk(6); ...
                          fbk(7)];
    ctrl_ref_curr= repmat(ctrl_ref,1,ctrl.N+1);
    % Solve the control problem
    [tau_l(i),tau_r(i),dyn_fbk_pred(:,i+1)]= controller(ctrl,sys,ctrl_ref_curr,fbk);

    fbk= propagate_plant(sys,fbk,[tau_l(i);tau_r(i)],ctrl.Ts);

end
main_loop_time = toc(main_loop);

ss_error = norm((fbk(1:3)-sim.xf(1:3)),2)

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
