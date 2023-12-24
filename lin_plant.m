function [state_next] = lin_plant(sys,state,ins,model)
    % This function implements the nonlinear plant model for the system.     

    state_lin = [state(4)+(sys.w/2)*state(5); ...
                 state(4)-(sys.w/2)*state(5); ...
                 state(6); ...
                 state(7)];

    state_lin_Dot= model.Ac*state_lin+model.Bc*ins;

    % State order X Y Psi v Psi_Dot Theta Theta_Dot
    v               = state(4);
    psi             = state(3);
    psi_dot         = state(5);
    theta_dot       = state_lin_Dot(3);
    v_dot           = (state_lin_Dot(1)+state_lin_Dot(2))/2;
    psi_ddot        = state_lin_Dot(2);
    theta_ddot      = state_lin_Dot(4);

    state_next= [v*cos(psi);v*sin(psi);psi_dot;v_dot;psi_ddot;theta_dot;theta_ddot];
end

