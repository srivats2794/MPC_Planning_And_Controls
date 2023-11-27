function [state_dot] = plant(sys,state,ins)
% This function implements the nonlinear plant model for the system.     

    tau_l= ins(1);
    tau_r= ins(2);
    
    % State order X Y Psi v Psi_Dot Theta Theta_Dot 
    v= state(4);
    psi= state(3);
    psi_dot= state(5);
    theta= state(6);
    theta_dot= state(7);

    Me= sys.m_b+2*sys.m_w+sys.m_p+(2*sys.j_w/(sys.r_w^2));
    Je = sys.j_m_th+sys.j_p_th;
    
    m_l_p= sys.m_p*sys.l;

    psi_term= ((sys.j_w+(sys.r_w^2)*sys.m_w)*sys.w^2)/(2*sys.r_w*sys.r_w);

    theta_ddot_num= ((m_l_p*cos(theta)/Me)* ...
        (tau_l/sys.r_w+tau_r/sys.r_w+m_l_p*(theta_dot^2)*sin(theta))) ...
        - m_l_p*sys.g*sin(theta);
    
    theta_ddot_den= ((m_l_p*cos(theta))^2/Me) - (sys.m_p*sys.l*sys.l+Je);

    v_dot_num= ((tau_l+tau_r)/sys.r_w) + m_l_p*(theta_dot^2)*sin(theta) - ...
        (((m_l_p^2)*sys.g*sin(theta)*cos(theta))/(sys.m_p*sys.l*sys.l+Je));

    v_dot_den= Me - (((m_l_p*cos(theta))^2)/(sys.m_p*sys.l*sys.l+Je));

    theta_ddot= theta_ddot_num/theta_ddot_den;

    v_dot= v_dot_num/v_dot_den;

    psi_ddot_const =  sys.w/(sys.j_psi+psi_term);
    
    psi_ddot_var= (tau_l-tau_r)/sys.r_w;
    
    psi_ddot= psi_ddot_const*psi_ddot_var;

    state_dot= [v*cos(psi);v*sin(psi);psi_dot;v_dot;psi_ddot;theta_dot;theta_ddot];
end

