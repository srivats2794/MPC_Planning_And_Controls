function [x0, u0] = shift_n(T, x0, u,f)
st = x0;
con = u(1,:)';
k1 = f(st, con);   % new 
    k2 = f(st + T/2*k1, con); % new
    k3 = f(st + T/2*k2, con); % new
    k4 = f(st + T*k3, con); % new
    st=st +T/6*(k1 +2*k2 +2*k3 +k4); % new    

x0 = full(st);

u0 = [u(2:size(u,1),:);u(size(u,1),:)];
end