function [state] = propagate_plant(sys,state,ins,Ts)
    % This function implements RK4 integrator to propagate plant dynamics
    % Used in combination with plant()

    k1 = plant(sys,state,ins);  
    k2 = plant(sys,state + Ts/2*k1, ins); 
    k3 = plant(sys,state + Ts/2*k2, ins); 
    k4 = plant(sys,state + Ts*k3, ins); 
    state=state + Ts/6*(k1 +2*k2 +2*k3 +k4); 

end

