function [grad_x, r, rd, rdd, ctrl_param, IC] = system3()

% x = [q1 q2 \dot{q1} \dot{q2}]

%% INITIAL CONDITION
% ********************************************************
IC.x = [1;-1;0;0];  
% IC.x = [1;1;0.3141;-0.3141];  
IC.u = [0; 0];

%% SYSTEM DYNAMICS
% ********************************************************
M = @(x) diag([20 30]);

C = @(x) diag([10 15]);
% C = @(x) (1+rand(2,2)) * 10;
    
G = @(x) [
    1; 1
    ];

d = @(x) zeros(2,1);
% d = @(x) [
%     tanh(x(1)) * 2
%     -tanh(x(2)) * 2
% ] * 1e2;

grad_x = @(x,u, t) blkdiag(eye(2), M(x)) \ ...
    ([zeros(2) eye(2); zeros(2) -C(x)] * x + [zeros(2);eye(2)] * (-G(x)+u+d(x)));

%% REFERENCE SIGNAL
% ********************************************************

% dt = 1e-4;
% r = @(x, t)[ % desired trajectory
%     heaviside(t - 3) + 1
%     heaviside(t - 3) - 1
%     ] ;
% 
% rd = @(x,t) [ 
%     1/dt * (heaviside(t - 3) - heaviside(t - 3-dt))
%     1/dt * (heaviside(t - 3) - heaviside(t - 3-dt))
% ];
% 
% rdd = @(x, t) [
%     0;0
% ];

r = @(x, t)[ % desired trajectory
    1+0.2*sin(0.5*pi*t)
    -1-0.2*sin(0.5*pi*t)
    ] ;

rd = @(x,t) [ 
    0.2*cos(0.5*pi*t)*0.5*pi
    -0.2*cos(0.5*pi*t)*0.5*pi
];

rdd = @(x, t) [
    -0.2*sin(0.5*pi*t)*(0.5*pi)^2
    +0.2*sin(0.5*pi*t)*(0.5*pi)^2
];

end