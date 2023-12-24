function [state] = propagate_plant(sys,state,ins,Ts,lin,model)
    % This function implements RK4 integrator to propagate plant dynamics
    % Used in combination with plant()
    if lin==0
        k1 = plant(sys,state,ins);
        k2 = plant(sys,state + Ts/2*k1, ins);
        k3 = plant(sys,state + Ts/2*k2, ins);
        k4 = plant(sys,state + Ts*k3, ins);
        state=state + Ts/6*(k1 +2*k2 +2*k3 +k4);
    else
        k1 = lin_plant(sys,state,ins,model);
        k2 = lin_plant(sys,state + Ts/2*k1, ins,model);
        k3 = lin_plant(sys,state + Ts/2*k2, ins,model);
        k4 = lin_plant(sys,state + Ts*k3, ins,model);
        state=state + Ts/6*(k1 +2*k2 +2*k3 +k4);
    end
end
