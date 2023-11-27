function [tau_l,tau_r,prediction] = controller(ctrl,sys,reference,feedback)
% This function updates the conic problem setup by OSQP and computes inputs

            %% Trasnforming measurement to controller states
            % Plant States -> X Y psi V psiDot theta thetaDot
            % Controller States -> X_L,X_R, theta, thetaDot
            x0= [feedback(4)+(sys.w/2)*feedback(5); ...
                 feedback(4)-(sys.w/2)*feedback(5); ...
                 feedback(6); ...
                 feedback(7)];
            
            %% State part of the objective function
            q_new = [reshape((ctrl.Q*reference),ctrl.N+1,1); zeros(ctrl.N*ctrl.nu, 1)];
            
            %% Equality constraints
            leq = [-x0; zeros(ctrl.N*ctrl.nx, 1)];
            ueq = leq;
            
            %% Inequality constraints -> input and state constraints
            lineq = [repmat(ctrl.x_min, ctrl.N+1, 1); repmat(ctrl.tau_min*ones(ctrl.nu,1), ctrl.N, 1)];
            uineq = [repmat(ctrl.x_max, ctrl.N+1, 1); repmat(ctrl.tau_max*ones(ctrl.nu,1), ctrl.N, 1)];
            % - OSQP constraints
            
            l_new = [leq; lineq];
            u_new = [ueq; uineq];

            ctrl.solver.update('q', q_new, 'l', l_new, 'u', u_new);
            
            %%%%%%% solve %%%%%%%
            res          = ctrl.solver.solve();

            % check solver status
            if ~strcmp(res.info.status, 'solved')
                error('OSQP did not solve the problem!')
            end

            %%%%%%% apply first control input to the plant %%%%%%%
            ins        = res.x((ctrl.N+1)*ctrl.nx+1:(ctrl.N+1)*ctrl.nx+ctrl.nu);
            tau_l= ins(1);
            tau_r= ins(2);
           
            % - MPC solution: x(k+1)
            prediction       = res.x(ctrl.nx+1:ctrl.nx+ctrl.nx);
end
