function [ctrl_sys] = ctrl_sys_setup(sys,ctrl)
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
    p3= (-p04*sys.g)/(Me*p01-p04);
    p4= p01/(Me*sys.r_w*p01-p04*sys.r_w);
    p5= (sys.w*sys.w*sys.r_w)/ ...
        ((2*sys.j_psi*sys.r_w*sys.r_w)+ ...
        ((sys.j_w+sys.r_w*sys.r_w*sys.m_w)*sys.w*sys.w));
    
    %% Order documentation
    % States order -> x_l_dot,x_r_dot,theta,theta_Dot
    % Inputs order -> tau_l,tau_r
    
    %% System matrices continuous time

    A= [0 0 p3 0 ;
        0 0 p3 0;
        0 0 1  0;
        0 0 p1 0];

    B= [p4+p5 p4-p5;
        p4-p5 p4+p5;
          0     0;
         p2    p2;];

    C= eye(4);

    D = zeros(4,2);
    
    %% Discretization

    sys= ss(A,B,C,D);
    sys_d= c2d(sys,ctrl.Ts);
    ctrl_sys.A= sys_d.A;
    ctrl_sys.B= sys_d.B;
    
end