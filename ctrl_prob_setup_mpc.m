function solver = ctrl_prob_setup_mpc(lin_sys,ctrl)
   
   P = blkdiag( kron(speye(ctrl.N), ctrl.Q), ctrl.Q, kron(speye(ctrl.N), ctrl.R) );
   q = [repmat(-ctrl.Q*zeros(4,1), ctrl.N+1,1); zeros(ctrl.N*ctrl.nu, 1)];
   Ax = kron(speye(ctrl.N+1), -speye(ctrl.nx)) + kron(sparse(diag(ones(ctrl.N, 1), -1)), lin_sys.Ad);
   Bu = kron([sparse(1, ctrl.N); speye(ctrl.N)], lin_sys.Bd);
   Aeq = [Ax, Bu];
   leq = [-zeros(4,1); zeros(ctrl.N*ctrl.nx, 1)];
   ueq = leq;
   % - input and state constraints
   Aineq = speye((ctrl.N+1)*ctrl.nx + ctrl.N*ctrl.nu);
   lineq = [repmat(ctrl.x_min, ctrl.N+1, 1); repmat(ctrl.tau_min*ones(ctrl.nu,1), ctrl.N, 1)];
   uineq = [repmat(ctrl.x_max, ctrl.N+1, 1); repmat(ctrl.tau_max*ones(ctrl.nu,1), ctrl.N, 1)];
   % - OSQP constraints
   A = [Aeq; Aineq];
   l = [leq; lineq];
   u = [ueq; uineq];

   % Create an OSQP object
   solver = osqp;

   % Setup workspace
   solver.setup(P, q, A, l, u, 'warm_start', true);

end

