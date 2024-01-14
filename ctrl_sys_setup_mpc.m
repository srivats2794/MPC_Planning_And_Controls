function [ctrl_sys] = ctrl_sys_setup_mpc(sys)
% This function computes the A and B matrix for the control design

    %% Predefining some repeated cluster terms
    Me= sys.m_b+2*sys.m_w+sys.m_p+((2*sys.j_w)/(sys.r_w^2));
    Je= sys.j_m_th+sys.j_p_th;
    
    %% Predifining some big constants to avoid typos
    p01= sys.m_p*sys.l*sys.l+Je;
    p02= sys.m_p*sys.g*sys.l;
    p03= sys.m_p*sys.l;
    p04= sys.m_p*sys.m_p*sys.l*sys.l;

    ctrl_sys.p1= (-p02*Me)/(p04-p01*Me);
    ctrl_sys.p2= p03/(sys.r_w*p04-Me*sys.r_w*p01);
    ctrl_sys.p3= (-p04*sys.g)/(Me*p01-p04);
    ctrl_sys.p4= p01/(Me*sys.r_w*p01-p04*sys.r_w);
    ctrl_sys.p5= (sys.w*sys.w*sys.r_w)/ ...
        ((2*sys.j_psi*sys.r_w*sys.r_w)+ ...
        ((sys.j_w+sys.r_w*sys.r_w*sys.m_w)*sys.w*sys.w));
    
    %% Order documentation
    % States order -> x_l_dot,x_r_dot,theta,theta_Dot
    % Inputs order -> tau_l,tau_r
    
    %% System matrices continuous time

    A= [0 0 ctrl_sys.p3 0 ;
        0 0 ctrl_sys.p3 0;
        0 0 0  1;
        0 0 ctrl_sys.p1 0];

    B= [ctrl_sys.p4+ctrl_sys.p5 ctrl_sys.p4-ctrl_sys.p5;
        ctrl_sys.p4-ctrl_sys.p5 ctrl_sys.p4+ctrl_sys.p5;
          0     0;
         ctrl_sys.p2    ctrl_sys.p2;];

    C= eye(4);

    D = zeros(4,2);
    
    %% Discretization

    ctrl_sys.A = A;
    ctrl_sys.B = B;
    ctrl_sys.C = C;
    ctrl_sys.D = D;
    
    ctrl_sys.A_th=     [0 1
                       ctrl_sys.p1 0];
    ctrl_sys.B_th=     [0 0;
                        ctrl_sys.p2 ctrl_sys.p2];
    ctrl_sys.C_th= eye(2);
    ctrl_sys.D_th= zeros(2,2);
    
end