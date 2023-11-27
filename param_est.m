function [xhat,pk] = param_est(xhat, pk, Nominal_parameters, EKF_parameters, v_kappa, feedbacks, un_KF)

%%%%%%%%%%%%%%%%%%%%%%%% Cornering Stiffness (CS) Estimator %%%%%%%%%%%%%%%%%%%%

    %******* This script provides a mechanism to estimate the cornering stiffnesses for use in the MPC of the lateral controller   
    %-- Input: 
    %         (1) estimated states at current step                              (xhat)
    %         (2) covariance matrix at current step                             (pk)
    %         (3) model parameters                                              (Nominal_parameters)
    %         (4) EKF parameters                                                (EKF_parameters)
    %         (5) reference velocity and curvature                              (v_kappa)
    %         (6) feedback from the truck                                       (feedbacks)
    %         (7) upper/lower bounds for CSs                                    (un_KF)
    %-- Output: 
    %         (1) estimated states at next step                                 (xhat)
    %         (2) covariance matrix at next step                                (pk)




%-- Standard deviation of process and measurement noises for EKF algorithm
% one gyroscope and one accelerometer
%        (1)first sensor measure directly the yaw rate r
%        (2) second one estimate the lateral acceleration

%-- uncertainty in the sensor measurements
% Sensor noise can be found in sensor specification sheet
% or it can be obtained from data analysis during steady-state manoeuvring

MNkf1                =EKF_parameters(1);   % for yaw rate (deg/s)
MNkf2                =EKF_parameters(2);   % for acceleration (m/s^2)

%-- plant uncertainty
% (1) steer angle sensor noise
% (2) cornering stiffness of the front and rear tires

SN1                 =EKF_parameters(3); % for steer angle (deg)
SN2                 =EKF_parameters(4); % for front cornering stiffness  (N/rad)
SN3                 =EKF_parameters(5);  % for rear cornering stiffness  (N/rad)

%-- sampling time
dt_control          =EKF_parameters(6);

%-- EKF
%--states & measurements & inputs
%     (i)   we have 6 states: [beta   r   betaDot  rDot  Cf     Cr]=[sleep_angle yaw_rate sleep_rate yaw_acceleration CS_f CS_r]
%     (ii)  we have 2 measurments: yaw rate and lateral acceleration
%     (iii) we have 1 input: tire angle


%-- other variables
Vx                 =v_kappa(1);
curr               =v_kappa(2);

%-- feedbacks
xDot_f            =feedbacks(1); % we use actual feedback velocity (xDot_f) instead of desired velocity (Vx)
yDot_f            =feedbacks(2);
phi_f             =feedbacks(3);
phiDot_f          =feedbacks(4);
X_f               =feedbacks(5);
Y_f               =feedbacks(6);
u                 =feedbacks(7); % tire angle

%-- model parameters
m_0           =Nominal_parameters(1);
l_f_0         =Nominal_parameters(2);
l_r_0         =Nominal_parameters(3);
I_z_0         =Nominal_parameters(4);
C_f_0         =Nominal_parameters(5);
C_r_0         =Nominal_parameters(6);


%-- number of tires
t_f          =2; % front
t_b          =4; % back

%-- Covariance of process and measurement noises for EKF algorithm
ynoisekf    =[MNkf1 MNkf2];
xdotnoisekf =[SN1 SN2 SN3];
Rkf         =diag([ynoisekf(1)^2 , ynoisekf(2)^2]);
Qkf         =diag([xdotnoisekf(1)^2 , xdotnoisekf(2)^2 , xdotnoisekf(3)^2]);

%-- measurement matrix
Hk         =[0 1 0 0 0 0;0 Vx Vx 0 0 0];


%% EKF algorithm
%-- Covariance matrix
Bwk      =[0                             0 0;...
    0                             0 0;...
    t_f*xhat(5)/(m_0*Vx)          0 0;...
    (t_f*xhat(5)*l_f_0)/I_z_0     0 0;...
    0                             1 0;...
    0                             0 1];

Qdk      =dt_control*Bwk*Qkf*Bwk';


%-- measurement output: yaw rate and lateral acceleration
ym         =[phiDot_f; yDot_f];

%-- Linearized system with nominal parameters for filter

A31       =-(xhat(5)+xhat(6))/(m_0*Vx);
A32       =-(((xhat(5)*l_f_0-xhat(6)*l_r_0)/(m_0*Vx^2))+1);
A35       =t_f*(+(u/(m_0*Vx)) - (xhat(1)/(m_0*Vx)) - ((l_f_0*xhat(2))/(m_0*Vx^2)));
A36       =t_b*(+((l_r_0*xhat(2))/(m_0*Vx^2)) -(xhat(1)/(m_0*Vx)));

A41       =-(xhat(5)*l_f_0 - xhat(6)*l_r_0)/(I_z_0);
A42       =-(xhat(5)*l_f_0^2 + xhat(6)*l_r_0^2)/(I_z_0*Vx);
A45       =t_f*(+((u*l_f_0)/I_z_0) - ((l_f_0^2*xhat(2))/(I_z_0*Vx)) - ((l_f_0*xhat(1))/I_z_0));
A46       =t_b*(+((l_r_0*xhat(1))/I_z_0) -((l_r_0^2*xhat(2))/(I_z_0*Vx)));

Ak        =[ 1,    0,  dt_control,     0,         0,   0;
    0,    1,       0,      dt_control,   0,   0;
    A31, A32,      0,         0,       A35, A36;
    A41, A42,      0,         0,       A45, A46;
    0,    0,       0,         0,         1,   0;
    0,    0,       0,         0,         0,   1];


%*** (I) priori estimation

%--estimated nonlinear system with nominal parameters

% xhat      =[beta_hat   r_hat    betaDot_hat   rDot_hat   Cf_hat      Cr_hat] \in R^6
% xhat      =[x_hat(1)  x_hat(2)   x_hat(3)     x_hat(4)   x_hat(5)    x_hat(6) ] \in R^6
xhat        =[xhat(1) + xhat(3)*dt_control;...
    xhat(2) + xhat(4)*dt_control;...
    -((t_f*xhat(5)+t_b*xhat(6))/(m_0*Vx))*xhat(1) - (((t_f*xhat(5)*l_f_0-t_b*xhat(6)*l_r_0)/(m_0*Vx^2))+1)*xhat(2) +  (t_f*xhat(5)*u)/(m_0*Vx);...
    -((t_f*xhat(5)*l_f_0-t_b*xhat(6)*l_r_0)/I_z_0)*xhat(1) - ((t_f*xhat(5)*l_f_0^2+t_b*xhat(6)*l_r_0^2)/(I_z_0*Vx))*xhat(2) +  (t_f*xhat(5)*u*l_f_0)/I_z_0;...
    xhat(5);...
    xhat(6)];

%-- discrete-time covariance
pk         =Ak*pk*Ak'+Qdk;

%*** (II) posteriori estimation

% --discrete kalman gain
k         =pk*Hk'/(Hk*pk*Hk'+Rkf);

% --discrete estimated state update
xhat      =xhat+k*(ym-Hk*xhat);

% discrete covariance
pk        =(eye(size(Qdk))-k*Hk)*pk;

%--Project state estimate onto inequality constraints for bounding the CSs

% --max/min values allowed
un_cf           = un_KF(1);
un_cr           = un_KF(2);
C_f_max         =C_f_0+C_f_0*un_cf;
C_f_min         =C_f_0-C_f_0*un_cf;
C_r_max         =C_r_0+C_r_0*un_cr;
C_r_min         =C_r_0-C_r_0*un_cr;

%--upper/lower bounds
lb              = [C_f_min; C_r_min];
ub              = [C_f_max; C_r_max];
Mx              = [0 0 0 0 1 0; 0 0 0 0 0 1;
    0 0 0 0 -1 0; 0 0 0 0 0 -1;];
dx              = [ub; -lb];

%-- projection operator
xhat             = proj(xhat, Mx, dx);

%-- post-opration for safety assurance
xhat(5)          = min(xhat(5),C_f_max);
xhat(5)          = max(xhat(5),C_f_min);
xhat(6)          = min(xhat(6),C_r_max);
xhat(6)          = max(xhat(6),C_r_min);

%% Projection operator

function xproj = proj(x, M, d)
% Project x onto the region defined by Mx <= d using the projected gradient method

xproj = x;

%-- define max iteration, step size, and tolerance
max_iter = 1e3;
step_size = 1;
tol = 1e-16;

%-- projection
for iter = 1:max_iter
    g = M'*(M*xproj - d);
    if norm(g, inf) < tol
        break;
    end
    xproj = xproj - step_size*g;
    xproj = min(max(xproj, min(x)), max(x));
end


