function [controller] = ctrl_setup_lqr(sys,ctrl)
% This function computes the A and B matrix for the control design

    %% Predefining some repeated cluster terms
    Me= sys.m_b+2*sys.m_w+sys.m_p+((2*sys.j_w)/(sys.r_w^2));
    Je= sys.j_m_th+sys.j_p_th;
    
    %% Predifining some big constants to avoid typos
    p01= sys.m_p*sys.l*sys.l+Je;
    p02= sys.m_p*sys.g*sys.l;
    p03= sys.m_p*sys.l;
    p04= sys.m_p*sys.m_p*sys.l*sys.l;

    p1= (-p02*Me)/(p04-p01*Me);
    p2= p03/(sys.r_w*p04-Me*sys.r_w*p01);
        
    %% Order documentation
    % States order -> theta,theta_Dot
    % Inputs order -> tau_l,tau_r
    
    %% System matrices continuous time

    A= [0  1;
        p1 0];

    B= [ 0     0;
         p2    p2;];

    C= eye(2);

    D = zeros(2,2);
    
    %% Discretization

    sys= ss(A,B,C,D);
    sys_d= c2d(sys,ctrl.Ts);
    controller.A= sys_d.A;
    controller.B= sys_d.B;

     [~,controller.K,~] = idare(controller.A,controller.B,ctrl.Q,ctrl.R,[],[]);

     controller.K= -controller.K;
end