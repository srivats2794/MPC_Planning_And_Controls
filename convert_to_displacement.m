function [x_d] = convert_to_displacement(X1,Y1,theta1,X0,Y0,theta0)
    e_th= theta1-theta0;

    rot = [cos(e_th) sin(e_th);
          -sin(e_th) cos(e_th)];

    body_coords= rot*[X1-X0;Y1-Y0];

    x= body_coords(1);
    y= body_coords(2);

    x_d= sqrt(x^2+y^2);
end

