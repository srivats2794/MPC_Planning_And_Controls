function [state_next] = lin_plant(sys,state,ins,model)
    % This function implements the nonlinear plant model for the system.     

  
    % This function implements the nonlinear plant model for the system.     

    v_l= state(4);
    v_r= state(5);    
    psi= state(3);
    theta= state(6);
    thetaDot= state(7);
    tau_L= ins(1);
    tau_R= ins(2);
     
    state_next= [ ((v_l+v_r)/2)*cos(psi); ...
                  ((v_l+v_r)/2)*sin(psi); ...
                  ((v_l-v_r)/sys.w); ...
                  (model.A*[v_l;v_r;theta;thetaDot]+model.B*[tau_L;tau_R])];

end

